#[cfg(not(target_arch = "wasm32"))]
pub fn available_workers() -> usize {
    std::thread::available_parallelism()
        .map(|value| value.get())
        .unwrap_or(1)
}

#[cfg(target_arch = "wasm32")]
pub fn available_workers() -> usize {
    1
}
