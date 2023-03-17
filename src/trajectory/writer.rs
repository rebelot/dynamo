use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;

use crate::system::System;

pub struct TrajectoryWriter {
    file: File,
    interval: i32,
}

impl TrajectoryWriter {
    pub fn new(filename: &str, interval: i32) -> Self {
        Self {
            file: OpenOptions::new()
                .append(true)
                .create(true)
                .open(filename)
                .unwrap(),
            interval,
        }
    }
    pub fn write(&mut self, system: &System, step: &i32, time: &f64) {
        if step % self.interval != 0 {
            return;
        }
        let mut pos = Vec::new();
        for a in system.topology.atoms.iter() {
            pos.push(a.pos[0].to_string());
            pos.push(a.pos[1].to_string());
            pos.push(a.pos[2].to_string());
        }

        if let Err(e) = writeln!(self.file, "{} {}", time, pos.join(" ")) {
            eprintln!("Couldn't write to file: {}", e);
        }
    }
}
