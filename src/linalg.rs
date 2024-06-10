use crate::{Rvec, DIM};

pub const DEG2RAD: f32 = std::f32::consts::PI / 180.0;

#[inline]
pub fn dot(avec: &Rvec, bvec: &Rvec) -> f32 {
    avec[0] * bvec[0] + avec[1] * bvec[1] + avec[2] * bvec[2]
    // return avec.iter().zip(bvec.iter()).map(|(a, b)| a * b).sum();
}

#[inline]
pub fn cross(avec: &Rvec, bvec: &Rvec) -> Rvec {
    [
        avec[1] * bvec[2] - avec[2] * bvec[1],
        avec[2] * bvec[0] - avec[0] * bvec[2],
        avec[0] * bvec[1] - avec[1] * bvec[0],
    ]
}

#[inline]
pub fn rmul(vec: &Rvec, x: f32) -> Rvec {
    let mut cvec = [0.0; DIM];
    for i in 0..DIM {
        cvec[i] = vec[i] * x;
    }
    cvec
}

#[inline]
pub fn rdiv(vec: &Rvec, x: f32) -> Rvec {
    let mut cvec = [0.0; DIM];
    for i in 0..DIM {
        cvec[i] = vec[i] / x;
    }
    cvec
}

#[inline]
pub fn radd(vec: &Rvec, x: f32) -> Rvec {
    let mut cvec = [0.0; DIM];
    for i in 0..DIM {
        cvec[i] = vec[i] + x;
    }
    cvec
}

#[inline]
pub fn rsub(vec: &Rvec, x: f32) -> Rvec {
    let mut cvec = [0.0; DIM];
    for i in 0..DIM {
        cvec[i] = vec[i] - x;
    }
    cvec
}

#[inline]
pub fn rvsub(avec: &Rvec, bvec: &Rvec) -> Rvec {
    let mut cvec = [0.0; DIM];
    for i in 0..DIM {
        cvec[i] = avec[i] - bvec[i];
    }
    cvec
}

#[inline]
pub fn rvadd(avec: &Rvec, bvec: &Rvec) -> Rvec {
    let mut cvec = [0.0; DIM];
    for i in 0..DIM {
        cvec[i] = avec[i] + bvec[i];
    }
    cvec
}

// pub fn angle(avec: &Rvec, bvec: &Rvec) -> f32 {
//     let x = dot(avec, bvec);
//     let y = norm2(&cross(avec, bvec)).sqrt();
//     return y.atan2(x);
// }

#[inline]
pub fn norm2(vec: &Rvec) -> f32 {
    dot(vec, vec)
}

#[inline]
pub fn min_image<'a>(vec: &'a mut Rvec, pbc: &Rvec) -> &'a Rvec {
    vec.iter_mut().zip(pbc.iter()).for_each(|(v, p)| {
        let half_box = 0.5 * p;
        if *v > half_box {
            *v -= p;
        } else if *v < -half_box {
            *v += p;
        }
    });
    vec
}

#[inline]
pub fn wrap<'a>(vec: &'a mut Rvec, pbc: &Rvec) -> &'a Rvec {
    vec.iter_mut().zip(pbc.iter()).for_each(|(v, p)| *v %= p);
    vec
}

#[inline]
pub fn displace_vec(avec: &Rvec, bvec: &Rvec) -> Rvec {
    let mut cvec = [0.0; DIM];
    // for i in 0..DIM {
    //     cvec[i] = bvec[i] - avec[i];
    // }
    cvec.iter_mut()
        .zip(avec.iter())
        .zip(bvec.iter())
        .for_each(|((c, a), b)| *c = b - a);
    cvec
}

// Iterator that return the indices (i, j) corresponding to elements
// of a lower-triangular matrix with shape n x n.
pub fn tril_indices_from(n: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..n).flat_map(move |i| (i..n).map(move |j| (i, j)))
}
