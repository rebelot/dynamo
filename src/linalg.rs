// use std::iter::Sum;
// use std::ops::Mul;

const DIM: usize = 3;

pub fn dot(avec: &[f64; DIM], bvec: &[f64; DIM]) -> f64 {
    // return avec[0] * bvec[0] + avec[1] * bvec[1] + avec[2] * bvec[2];
    return avec.iter().zip(bvec.iter()).map(|(a, b)| a * b).sum();
}

// pub fn dot<T>(avec: &Vec<T>, bvec: &Vec<T>) -> T
// where
//     T: Sum + Mul<T, Output = T> + Copy,
// {
//     return avec.iter().zip(bvec.iter()).map(|(a, b)| *a * *b).sum();
// }

// pub fn cos_angle(avec: &[f64; DIM], bvec: &[f64; DIM]) -> f64 {
//     let ab = dot(avec, bvec);
//     let anorm = dot(avec, avec).sqrt();
//     let bnorm = dot(bvec, bvec).sqrt();
//     return ab / (anorm * bnorm);
// }

pub fn cross(avec: &[f64; 3], bvec: &[f64; 3]) -> [f64; 3] {
    return [
        avec[1] * bvec[2] - avec[2] * bvec[1],
        avec[2] * bvec[0] - avec[0] * bvec[2],
        avec[0] * bvec[1] - avec[1] * bvec[0],
    ];
}

pub fn angle(avec: &[f64; DIM], bvec: &[f64; DIM]) -> f64 {
    let x = dot(avec, bvec);
    let y = norm2(&cross(avec, bvec)).sqrt();
    return y.atan2(x);
}

pub fn norm2(vec: &[f64; DIM]) -> f64 {
    return dot(vec, vec);
}

pub fn min_image<'a>(vec: &'a mut [f64; DIM], pbc: &[f64; DIM]) -> &'a [f64; DIM] {
    for i in 0..DIM {
        let half_box = 0.5 * pbc[i];
        if vec[i] > half_box {
            vec[i] -= pbc[i];
        } else if vec[i] < -half_box {
            vec[i] += pbc[i];
        }
    }
    return vec;
}

pub fn wrap<'a>(vec: &'a mut [f64; DIM], pbc: &[f64; DIM]) -> &'a [f64; DIM] {
    for i in 0..DIM {
        vec[i] %= pbc[i];
    }
    return vec;
}

pub fn displace_vec(avec: &[f64; DIM], bvec: &[f64; DIM]) -> [f64; DIM] {
    let mut cvec = [0.0; DIM];
    for i in 0..DIM {
        cvec[i] = bvec[i] - avec[i];
    }
    return cvec;
}
