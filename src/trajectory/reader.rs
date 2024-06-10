use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::Rvec;

pub fn read_coords(file: &str) -> (Rvec, Vec<Rvec>) {
    let mut coords = Vec::new();
    let file = File::open(file).unwrap();
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        let mut xyz = [0.0; 3];
        line.split_whitespace()
            .enumerate()
            .for_each(|(i, x)| xyz[i] = x.parse::<f32>().unwrap());
        coords.push([xyz[0], xyz[1], xyz[2]])
    }
    (coords[0], coords[1..].to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_coords() {
        let (ref pbc, ref coords) = read_coords("tests/diala.crd");
        println!("{:?}", coords);
    }
}
