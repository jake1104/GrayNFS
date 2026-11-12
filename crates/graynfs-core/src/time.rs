pub use std::time::Duration;
#[cfg(not(target_arch = "wasm32"))]
pub use std::time::Instant;

#[cfg(target_arch = "wasm32")]
pub use wasm::Instant;

#[cfg(target_arch = "wasm32")]
mod wasm {
    use std::cmp::Ordering;
    use std::ops::{Add, Sub};
    use std::time::Duration;

    #[derive(Clone, Copy, Debug)]
    pub struct Instant(f64);

    fn monotonic_ms() -> f64 {
        use wasm_bindgen::{JsCast, JsValue};
        let global = js_sys::global();
        if let Ok(performance) = js_sys::Reflect::get(&global, &JsValue::from_str("performance"))
            && let Ok(now) = js_sys::Reflect::get(&performance, &JsValue::from_str("now"))
            && let Some(function) = now.dyn_ref::<js_sys::Function>()
            && let Ok(value) = function.call0(&performance)
            && let Some(ms) = value.as_f64()
        {
            return ms;
        }
        js_sys::Date::now()
    }

    impl Instant {
        pub fn now() -> Self {
            Self(monotonic_ms())
        }

        pub fn elapsed(&self) -> Duration {
            Instant::now().saturating_duration_since(*self)
        }

        pub fn duration_since(&self, earlier: Instant) -> Duration {
            Duration::from_secs_f64(((self.0 - earlier.0) / 1000.0).max(0.0))
        }

        pub fn saturating_duration_since(&self, earlier: Instant) -> Duration {
            self.duration_since(earlier)
        }

        pub fn checked_add(&self, duration: Duration) -> Option<Self> {
            Some(Self(self.0 + duration.as_secs_f64() * 1000.0))
        }

        pub fn checked_sub(&self, duration: Duration) -> Option<Self> {
            Some(Self(self.0 - duration.as_secs_f64() * 1000.0))
        }
    }

    impl PartialEq for Instant {
        fn eq(&self, other: &Self) -> bool {
            self.0.total_cmp(&other.0) == Ordering::Equal
        }
    }

    impl Eq for Instant {}

    impl PartialOrd for Instant {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for Instant {
        fn cmp(&self, other: &Self) -> Ordering {
            self.0.total_cmp(&other.0)
        }
    }

    impl Add<Duration> for Instant {
        type Output = Instant;

        fn add(self, duration: Duration) -> Instant {
            Self(self.0 + duration.as_secs_f64() * 1000.0)
        }
    }

    impl Sub<Duration> for Instant {
        type Output = Instant;

        fn sub(self, duration: Duration) -> Instant {
            Self(self.0 - duration.as_secs_f64() * 1000.0)
        }
    }

    impl Sub<Instant> for Instant {
        type Output = Duration;

        fn sub(self, earlier: Instant) -> Duration {
            self.duration_since(earlier)
        }
    }

    #[cfg(test)]
    mod tests {
        use super::Instant;
        use std::time::Duration;
        use wasm_bindgen_test::wasm_bindgen_test;

        #[wasm_bindgen_test]
        fn elapsed_advances_above_zero() {
            let start = Instant::now();
            let mut spins: u32 = 0;
            while start.elapsed().is_zero() {
                spins += 1;
                assert!(spins < 10_000_000);
            }
        }

        #[wasm_bindgen_test]
        fn backend_is_monotonic_tick_not_wall_clock() {
            let tick = Instant::now();
            assert!(tick.0 < js_sys::Date::now());
            assert!(tick.0 >= 0.0);
        }

        #[wasm_bindgen_test]
        fn now_is_monotonic_non_decreasing() {
            let first = Instant::now();
            let second = Instant::now();
            assert!(second >= first);
            assert!(second.saturating_duration_since(first) < Duration::from_secs(60));
            assert!(first.saturating_duration_since(second).is_zero());
        }

        #[wasm_bindgen_test]
        fn add_sub_round_trip() {
            let start = Instant::now();
            let delta = Duration::from_millis(1500);
            let there = start + delta;
            assert!(there > start);
            assert!(start.checked_add(delta).is_some());
            assert!(there.checked_sub(delta).is_some());
            let back = there - delta;
            assert!(back < there);
            let error_ms = if back.0 >= start.0 {
                back.0 - start.0
            } else {
                start.0 - back.0
            };
            assert!(error_ms < 0.001);
            let measured = there - start;
            assert!(measured <= delta + Duration::from_micros(10));
            assert!(measured + Duration::from_micros(10) >= delta);
        }
    }
}
