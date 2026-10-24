# Local gate mirroring .github/workflows/ci.yml step order.
# Variants: plain cargo test (nextest absent); full-workspace release build.
param(
  [switch]$Full,
  [string]$LogRoot = "logs/gate"
)

$ErrorActionPreference = "Continue"

$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $RepoRoot

$stamp = Get-Date -Format "yyyy.MM.dd-HHmmss"
$gateDir = Join-Path $LogRoot $stamp
New-Item -ItemType Directory -Force -Path $gateDir | Out-Null
$summaryFile = Join-Path $gateDir "summary.log"

$probe = Get-Command nextest -ErrorAction SilentlyContinue
if ($probe -eq $null) {
  $nextestState = "absent"
} else {
  $nextestState = "present"
}

"gate start {0} Full={1}" -f $stamp, $Full | Tee-Object -FilePath $summaryFile
"variant nextest={0} so test step uses plain cargo test --workspace" -f $nextestState | Add-Content -Path $summaryFile
"variant release-build uses full-workspace cargo build --release (CI uses -p graynfs-cli)" | Add-Content -Path $summaryFile

function Invoke-GateStep {
  param(
    [string]$StepName,
    [scriptblock]$Action
  )
  $logFile = Join-Path $gateDir ($StepName + ".log")
  New-Item -ItemType File -Force -Path $logFile | Out-Null
  Write-Host ("[gate] step {0}" -f $StepName)
  $started = Get-Date
  & $Action 2>&1 | Tee-Object -FilePath $logFile
  $code = $LASTEXITCODE
  $secs = [math]::Round(((Get-Date) - $started).TotalSeconds, 1)
  $line = "{0} exit={1} elapsed={2}s" -f $StepName, $code, $secs
  $line | Add-Content -Path $summaryFile
  Write-Host ("[gate] {0}" -f $line)
  if ($code -ne 0) {
    Write-Host ("[gate] FAIL at {0} code {1}" -f $StepName, $code)
    exit $code
  }
}

Invoke-GateStep "fmt" { cargo fmt --check }
Invoke-GateStep "clippy" { cargo clippy --workspace --all-targets -- -D warnings }
Invoke-GateStep "test" { cargo test --workspace }
Invoke-GateStep "release-build" { cargo build --release }

if ($Full) {
  # quarantine-105bit policy mirror (owner CI-quarantine, review 2026-10-20):
  # blocking-transition gate v0.1.4.3.2: B1 5 consecutive release greens inside
  # 5-20s band (v0.1.4.3.1 7.48s = point 1/5) + B2 one debug green within 1800s
  # byte-identical factors. Until B1+B2 met, DECLINE (CI continue-on-error stays
  # true; this local step/filter unchanged). Filter append-only, never narrow.
  # verdict v0.1.4.5.3 (review 2026-10-20): B1+B2 BOTH MET (see ci.yml verdict
  # block for counts). Local step/filter unchanged here; non-blocking stays
  # pending the separate gated blocking-flip decision (READY-TO-TRANSITION).
  # verdict v0.1.4.8.1 (review 2026-10-20): blocking flip DECLINED (default);
  # READY-TO-TRANSITION maintained, NO flip (see ci.yml verdict block).
  Invoke-GateStep "quarantine-105bit" { cargo test --release -p graynfs-core --lib factors_105bit -- --ignored }
  Invoke-GateStep "wheel-build" { maturin build --release --out wheels-gate }
  $wheel = Get-ChildItem -Path "wheels-gate" -Filter "graynfs-*.whl" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
  if ($wheel -eq $null) {
    "wheel-install exit=1 elapsed=0s (no wheel found)" | Add-Content -Path $summaryFile
    Write-Host "[gate] FAIL at wheel-install code 1 (no wheel found)"
    exit 1
  }
  $wheelPath = $wheel.FullName
  Invoke-GateStep "wheel-install" { pip install --force-reinstall --no-deps $wheelPath }
  Invoke-GateStep "pytest" { python -m pytest -q python/tests }
  Invoke-GateStep "version-guard" { python scripts/check_version.py }
  # quintic-scope (2026.10.19.7.3): quintic pin is 27-bit 100160063 only;
  # 43-bit row 4606404410723 deferred (efficiency-only, 60 min budget).
  # wide_gnfs expects siqs; filter append-only, never narrow.
  $degreesProbe = Get-Command cargo-nextest -ErrorAction SilentlyContinue
  if ($degreesProbe -ne $null) {
    Invoke-GateStep "gnfs-degrees" { cargo nextest run --release -p graynfs-core --run-ignored=only -E 'test(cubic_gnfs) + test(quartic_gnfs) + test(quintic_gnfs) + test(wide_gnfs)' }
  } else {
    Invoke-GateStep "gnfs-degrees" { cargo test --release -p graynfs-core --test core -- --ignored cubic_gnfs quartic_gnfs quintic_gnfs wide_gnfs }
  }
  Invoke-GateStep "perf-smoke" { python scripts/perf_smoke.py }
}

"gate PASS" | Add-Content -Path $summaryFile
Write-Host "[gate] PASS"
exit 0
