use perfcnt::{
    AbstractPerfCounter, PerfCounter,
    linux::{PerfCounterBuilderLinux, SoftwareEventType},
};
use std::io;
use std::time::{Duration, Instant};

pub struct StageTiming {
    start: Instant,
    task_clocks: Option<Vec<PerfCounter>>,
}

impl StageTiming {
    pub fn start() -> Self {
        // A single inherited counter misses threads in GGCAT's persistent Rayon pool.
        // Attach to every existing thread, and inherit counters for threads they spawn.
        let task_clocks = (|| -> io::Result<Vec<PerfCounter>> {
            let tids = std::fs::read_dir("/proc/self/task")?
                .map(|entry| {
                    entry?
                        .file_name()
                        .to_string_lossy()
                        .parse::<i32>()
                        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
                })
                .collect::<io::Result<Vec<_>>>()?;
            tids.into_iter()
                .map(|tid| {
                    let counter =
                        PerfCounterBuilderLinux::from_software_event(SoftwareEventType::TaskClock)
                            .for_pid(tid)
                            .inherit()
                            .disable()
                            .finish()?;
                    counter.start()?;
                    Ok(counter)
                })
                .collect()
        })()
        .map_err(|error| eprintln!("CPU time measurement unavailable: {error}"))
        .ok();
        Self {
            start: Instant::now(),
            task_clocks,
        }
    }

    pub fn finish(mut self) -> String {
        let wall = self.start.elapsed();
        let cpu = self.task_clocks.as_mut().and_then(|counters| {
            counters
                .iter_mut()
                .map(|counter| {
                    counter.stop()?;
                    counter.read()
                })
                .sum::<io::Result<u64>>()
                .map(Duration::from_nanos)
                .map_err(|error| eprintln!("CPU time measurement unavailable: {error}"))
                .ok()
        });
        match cpu {
            Some(cpu) => format!(
                "{wall:.2?} wall, {cpu:.2?} CPU, {:.1}% average CPU",
                100.0 * cpu.as_secs_f64() / wall.as_secs_f64()
            ),
            None => format!("{wall:.2?} wall, CPU unavailable"),
        }
    }
}
