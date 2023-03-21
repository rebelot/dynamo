use crate::ffield::*;

use super::{System, Topology};
use std::fs::File;
use std::io::{BufRead, BufReader};

enum Section {
    MAIN,
    NONBOND,
    BOND,
    ANGLE,
    DIHEDRAL,
    BOX,
}

impl Default for Section {
    fn default() -> Self {
        Section::MAIN
    }
}

pub fn read(filename: &str) -> System {
    let mut top = Topology::new();
    let mut forces = Forces::new();
    let f = File::open(filename).unwrap();

    let mut section: Section = Default::default();
    let mut section_counter = 0;

    let pbc = &mut [0.0, 0.0, 0.0];

    for line in BufReader::new(f).lines() {
        if let Ok(line) = line {
            if line.is_empty() {
                continue;
            }
            if line.starts_with("MAIN") {
                section = Section::MAIN;
                section_counter = 0;
            } else if line.starts_with("NONBOND") {
                section = Section::NONBOND;
                section_counter = 0;
            } else if line.starts_with("BOND") {
                section = Section::BOND;
                section_counter = 0;
            } else if line.starts_with("ANGLE") {
                section = Section::ANGLE;
                section_counter = 0;
            } else if line.starts_with("DIHEDRAL") {
                section = Section::DIHEDRAL;
                section_counter = 0;
            } else if line.starts_with("BOX") {
                section = Section::BOX;
                section_counter = 0;
            } else {
                section_counter += 1;
            }

            if section_counter > 0 {
                let fields: Vec<&str> = line.split_whitespace().collect();

                match section {
                    Section::MAIN => {
                        let atomtype = fields[0].to_string();
                        let name = fields[1].to_string();
                        let element = fields[2].to_string();
                        let mass = fields[3].parse::<f32>().unwrap();
                        let vdw = fields[4].parse::<f32>().unwrap();
                        let charge = fields[5].parse::<f32>().unwrap();
                        let x = fields[6].parse::<f32>().unwrap();
                        let y = fields[7].parse::<f32>().unwrap();
                        let z = fields[8].parse::<f32>().unwrap();
                        top.add_atom(atomtype, name, element, mass, vdw, charge, [x, y, z]);
                    }
                    Section::NONBOND => {
                        let a = fields[0].parse::<usize>().unwrap();
                        let c12 = fields[1].parse::<f32>().unwrap();
                        let c6 = fields[2].parse::<f32>().unwrap();
                        // top.atoms[a].lj.c12 = c12;
                        // top.atoms[a].lj.c6 = c6;
                    }
                    Section::BOND => {
                        let a1 = fields[0].parse::<usize>().unwrap();
                        let a2 = fields[1].parse::<usize>().unwrap();
                        let k = fields[2].parse::<f32>().unwrap();
                        let r0 = fields[3].parse::<f32>().unwrap();
                        forces.bonds.push(BondHarmonic::new(k, r0, [a1, a2]))
                    }
                    Section::ANGLE => {
                        let a1 = fields[0].parse::<usize>().unwrap();
                        let a2 = fields[1].parse::<usize>().unwrap();
                        let a3 = fields[2].parse::<usize>().unwrap();
                        let k = fields[3].parse::<f32>().unwrap();
                        let t0 = fields[4].parse::<f32>().unwrap();
                        forces.angles.push(AngleHarmonic::new(k, t0, [a1, a2, a3]))
                    }
                    Section::DIHEDRAL => {
                        let a1 = fields[0].parse::<usize>().unwrap();
                        let a2 = fields[1].parse::<usize>().unwrap();
                        let a3 = fields[2].parse::<usize>().unwrap();
                        let a4 = fields[3].parse::<usize>().unwrap();
                        let k = fields[4].parse::<f32>().unwrap();
                        let p0 = fields[5].parse::<f32>().unwrap();
                        let n = fields[6].parse::<f32>().unwrap();
                        forces
                            .torsions
                            .push(DihedralPeriodic::new(k, p0, n, [a1, a2, a3, a4]))
                    }
                    Section::BOX => {
                        pbc[0] = fields[0].parse::<f32>().unwrap();
                        pbc[1] = fields[1].parse::<f32>().unwrap();
                        pbc[1] = fields[2].parse::<f32>().unwrap();
                    }
                }
            }
        }
    }
    return System::new(top, forces, *pbc);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_reads() {
        let sys = read("tests/simple.sys");
    }
}
