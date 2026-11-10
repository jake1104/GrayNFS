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

    impl Instant {
        pub fn now() -> Self {
            Self(js_sys::Date::now())
        }

        pub fn elapsed(&self) -> Duration {
            self.saturating_duration_since(*self)
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
}
