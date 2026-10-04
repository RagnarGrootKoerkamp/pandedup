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

/// Peak resident memory of the process. On Linux, `ru_maxrss` is in KiB;
/// on macOS, it is in bytes.
fn process_max_rss() -> io::Result<u64> {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
        return Err(io::Error::last_os_error());
    }
    let usage = unsafe { usage.assume_init() };
    let rss = u64::try_from(usage.ru_maxrss)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "negative ru_maxrss"))?;
    Ok(if cfg!(target_os = "macos") {
        rss
    } else {
        rss.saturating_mul(1024)
    })
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
        let max_rss = process_max_rss()
            .map_err(|error| warn!("MaxRSS measurement unavailable: {error}"))
            .ok();
        let mut details = Vec::with_capacity(2);
        if let Some(cpu) = cpu {
            details.push(format!(
                "{cpu:.2?} CPU; {:.1}%",
                100.0 * cpu.as_secs_f64() / wall.as_secs_f64()
            ));
        }
        if let Some(max_rss) = max_rss {
            details.push(if max_rss >= 1 << 30 {
                format!("MaxRSS {:.2} GiB", max_rss as f64 / (1u64 << 30) as f64)
            } else {
                format!("MaxRSS {:.1} MiB", max_rss as f64 / (1u64 << 20) as f64)
            });
        }
        if details.is_empty() {
            format!("{wall:.2?} wall")
        } else {
            format!("{wall:.2?} wall [{}]", details.join("; "))
        }
    }
}
