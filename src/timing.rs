use std::io;
use std::time::{Duration, Instant};
use tracing::warn;

fn process_cpu_time() -> io::Result<Duration> {
    let mut time = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    // CLOCK_PROCESS_CPUTIME_ID includes CPU used by every thread in this process.
    if unsafe { libc::clock_gettime(libc::CLOCK_PROCESS_CPUTIME_ID, &mut time) } != 0 {
        return Err(io::Error::last_os_error());
    }
    Ok(Duration::new(time.tv_sec as u64, time.tv_nsec as u32))
}

pub struct StageTiming {
    start: Instant,
    cpu_start: Option<Duration>,
}

impl StageTiming {
    pub fn start() -> Self {
        let cpu_start = process_cpu_time()
            .map_err(|error| warn!("CPU time measurement unavailable: {error}"))
            .ok();
        Self {
            start: Instant::now(),
            cpu_start,
        }
    }

    pub fn finish(self) -> String {
        let wall = self.start.elapsed();
        let cpu = self.cpu_start.and_then(|start| {
            process_cpu_time()
                .map(|end| end.saturating_sub(start))
                .map_err(|error| warn!("CPU time measurement unavailable: {error}"))
                .ok()
        });
        match cpu {
            Some(cpu) => format!(
                "{wall:.2?} wall [{cpu:.2?} CPU; {:.1}%]",
                100.0 * cpu.as_secs_f64() / wall.as_secs_f64()
            ),
            None => format!("{wall:.2?} wall"),
        }
    }
}
