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
  $degreesProbe = Get-Command cargo-nextest -ErrorAction SilentlyContinue
  if ($degreesProbe -ne $null) {
    Invoke-GateStep "gnfs-degrees" { cargo nextest run --release -p graynfs-core --run-ignored=only -E 'test(cubic_gnfs) + test(quartic_gnfs) + test(wide_gnfs)' }
  } else {
    Invoke-GateStep "gnfs-degrees" { cargo test --release -p graynfs-core --test core -- --ignored cubic_gnfs quartic_gnfs wide_gnfs }
  }
}

"gate PASS" | Add-Content -Path $summaryFile
Write-Host "[gate] PASS"
exit 0
