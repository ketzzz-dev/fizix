mod world;
mod body;
mod constraint;

pub type Precision = f64;

pub const EPSILON: Precision = 1e-9;
pub const EPSILON_SQUARED: Precision = EPSILON * EPSILON;

pub fn get_two_mut<T>(arr: &mut [T], i: usize, j: usize) -> Option<(&mut T, &mut T)> {
    if i == j || i >= arr.len() || j >= arr.len() {
        return None;
    }

    let (a, b) = if i < j {
        let (head, tail) = arr.split_at_mut(j);

        (&mut head[i], &mut tail[0])
    } else {
        let (head, tail) = arr.split_at_mut(i);
        
        (&mut tail[0], &mut head[j])
    };

    Some((a, b))
}

pub use world::*;
pub use body::*;
pub use constraint::*;