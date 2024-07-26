use ark_ff::{batch_inversion, Field, Fp64, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "18446744069414584321"]
#[generator = "7"]
pub struct FrConfig64;
pub type Field64 = Fp64<MontBackend<FrConfig64, 1>>;

pub fn eq_poly<F>(coords: &[F], mut point: usize) -> F
where
    F: Field,
{
    let n_variables = coords.len();
    assert!(point < (1 << n_variables));

    let mut acc = F::ONE;

    for val in coords.iter().rev() {
        let b = point % 2;
        acc *= if b == 1 { *val } else { F::ONE - *val };
        point >>= 1;
    }

    acc
}

pub struct LagrangePolynomialGray<F: Field> {
    position_bin: usize,
    position_gray: usize,
    value: F,
    precomputed: Vec<F>,
    num_variables: usize,
}

impl<F: Field> LagrangePolynomialGray<F> {
    pub fn new(point: &[F]) -> Self {
        let num_variables = point.len();
        assert!(point.iter().all(|&p| p != F::ZERO && p != F::ONE));

        // This is negated[i] = eq_poly(z_i, 0) = 1 - z_i
        let negated_points: Vec<_> = point.iter().map(|z| F::ONE - z).collect();
        // This is points[i] = eq_poly(z_i, 1) = z_i
        let points = point.to_vec();

        let mut to_invert = [negated_points.clone(), points.clone()].concat();
        batch_inversion(&mut to_invert);

        let (denom_0, denom_1) = to_invert.split_at(num_variables);

        let mut precomputed = vec![F::ZERO; 2 * num_variables];
        for n in 0..num_variables {
            precomputed[2 * n] = points[n] * denom_0[n];
            precomputed[2 * n + 1] = negated_points[n] * denom_1[n];
        }

        LagrangePolynomialGray {
            position_gray: gray_encode(0),
            position_bin: 0,
            value: negated_points.into_iter().product(),
            num_variables,
            precomputed,
        }
    }
}

pub fn gray_encode(integer: usize) -> usize {
    (integer >> 1) ^ integer
}

pub fn gray_decode(integer: usize) -> usize {
    match integer {
        0 => 0,
        _ => integer ^ gray_decode(integer >> 1),
    }
}

impl<F: Field> Iterator for LagrangePolynomialGray<F> {
    type Item = (usize, F);

    fn next(&mut self) -> Option<Self::Item> {
        if self.position_bin >= (1 << self.num_variables) {
            return None;
        }

        let result = (self.position_gray, self.value);

        let prev = self.position_gray;

        self.position_bin += 1;
        self.position_gray = gray_encode(self.position_bin);

        if self.position_bin < (1 << self.num_variables) {
            let diff = prev ^ self.position_gray;
            let i = (self.num_variables - 1) - diff.trailing_zeros() as usize;
            let flip = (diff & self.position_gray == 0) as usize;

            self.value *= self.precomputed[2 * i + flip];
        }

        Some(result)
    }
}

fn main() {}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use crate::Field64;

    use super::*;

    type F = Field64;

    #[test]
    fn test_gray_ordering() {
        let values = [
            (0b0000, 0b0000),
            (0b0001, 0b0001),
            (0b0010, 0b0011),
            (0b0011, 0b0010),
            (0b0100, 0b0110),
            (0b0101, 0b0111),
            (0b0110, 0b0101),
            (0b0111, 0b0100),
            (0b1000, 0b1100),
            (0b1001, 0b1101),
            (0b1010, 0b1111),
            (0b1011, 0b1110),
            (0b1100, 0b1010),
            (0b1101, 0b1011),
            (0b1110, 0b1001),
            (0b1111, 0b1000),
        ];

        for (bin, gray) in values {
            assert_eq!(gray_encode(bin), gray);
            assert_eq!(gray_decode(gray), bin);
        }
    }

    #[test]
    fn test_gray_ordering_iterator() {
        let point = vec![F::from(2), F::from(3), F::from(4)];

        for (i, (b, _)) in LagrangePolynomialGray::new(&point).enumerate() {
            assert_eq!(b, gray_encode(i));
        }
    }

    #[test]
    fn test_gray() {
        let point = vec![F::from(2), F::from(3), F::from(4)];

        let eq_poly_res: BTreeSet<_> = (0..(1 << 3)).map(|b| (b, eq_poly(&point, b))).collect();

        let gray_res: BTreeSet<_> = LagrangePolynomialGray::new(&point).collect();

        assert_eq!(eq_poly_res, gray_res);
    }

    #[test]
    fn test_gray_again() {
        let point = vec![F::from(2), F::from(3), F::from(4), F::from(5)];

        let eq_poly_res: BTreeSet<_> = (0..(1 << 4)).map(|b| (b, eq_poly(&point, b))).collect();

        let gray_res: BTreeSet<_> = LagrangePolynomialGray::new(&point).collect();

        assert_eq!(eq_poly_res, gray_res);
    }
}
