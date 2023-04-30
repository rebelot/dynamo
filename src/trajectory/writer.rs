use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;

use crate::system::System;
use crate::Rvec;

pub struct TrajectoryWriter {
    file: File,
    interval: i32,
}

impl TrajectoryWriter {
    pub fn new(filename: &str, interval: i32, append: bool) -> Self {
        Self {
            file: OpenOptions::new()
                .create(true)
                .write(true)
                .append(append)
                .open(filename)
                .unwrap(),
            interval,
        }
    }
    pub fn write(&mut self, positions: &[Rvec], step: &i32, time: &f32) {
        if step % self.interval != 0 {
            return;
        }

        let blob = positions
            .iter()
            .flat_map(|xyz| xyz.map(|p| p.to_string()))
            .collect::<Vec<String>>()
            .join(" ");

        if let Err(e) = writeln!(self.file, "{} {}", time, blob) {
            eprintln!("Couldn't write to file: {}", e);
        }
        self.file.flush().unwrap();
    }
}
