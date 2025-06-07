/// This module provides functionalities related to seed-based random number generation.
///
/// It includes:
/// - `events`: Defines events related to seed changes.
/// - `plugin`: Provides a Bevy plugin for integrating the seed resource and handling seed events.
/// - `resource`: Defines the `SeedResource` which manages the seed and random number generation logic.
///
/// It also exports several macros for convenient random number generation:
/// - `next_random_value!`: Generates the next random `f32` value, optionally for a specific category.
/// - `next_random_normalized!`: Generates the next random `f32` value normalized to the range `[0.0, 1.0)`, optionally for a specific category.
/// - `next_random_bool!`: Generates the next random `bool` value (approximately 50/50 chance), optionally for a specific category.
/// - `next_random_ratio!`: Generates the next random `bool` based on a given probability, optionally for a specific category. (Note: The macro name might be misleading, consider `next_random_bool_with_probability!` for clarity).
/// - `next_random_iter!`: Generates the next random `i64` within a specified range (exclusive upper bound), optionally for a specific category.
/// - `next_random_range!`: Generates the next random `f32` within a specified range (exclusive upper bound), optionally for a specific category.
pub mod events;
pub mod plugin;
pub mod resource;

/// Generates the next random `f32` value using the `SeedResource`.
///
/// The noise values are typically in the range `[-1.0, 1.0]`.
///
/// # Usage
///
/// ## With a specific category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`
/// // and `MyCategory` is a string or a type that can be converted to a string.
/// // let random_val = next_random_value!(seed_resource, "MyCategory");
/// ```
///
/// ## With the default category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`
/// // let random_val = next_random_value!(seed_resource);
/// ```
#[macro_export]
macro_rules! next_random_value {
    // Variante für spezifische Kategorie
    ($seed_resource_mut:expr, $category:expr) => {
        $seed_resource_mut.next_value_for_category($category)
    };
    // Variante für Standard/Default-Wert (muss NACH der spezifischeren Variante stehen)
    ($seed_resource_mut:expr) => {
        $seed_resource_mut.next_value()
    };
}

/// Generates the next random `f32` value, normalized to the range `[0.0, 1.0)`, using the `SeedResource`.
///
/// # Usage
///
/// ## With a specific category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`
/// // and `MyCategory` is a string or a type that can be converted to a string.
/// // let normalized_val = next_random_normalized!(seed_resource, "MyCategory");
/// ```
///
/// ## With the default category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`
/// // let normalized_val = next_random_normalized!(seed_resource);
/// ```
#[macro_export]
macro_rules! next_random_normalized {
    // Variante für spezifische Kategorie
    ($seed_resource_mut:expr, $category:expr) => {
        $seed_resource_mut.next_value_normalized(Some($category))
    };
    // Variante für Standard/Default-Wert (muss NACH der spezifischeren Variante stehen)
    ($seed_resource_mut:expr) => {
        $seed_resource_mut.next_value_normalized(None)
    };
}

/// Generates the next random `bool` value (approximately 50/50 chance) using the `SeedResource`.
///
/// This is achieved by checking if the underlying noise value (typically in `[-1.0, 1.0]`) is greater than `0.0`.
///
/// # Usage
///
/// ## With a specific category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`
/// // and `MyCategory` is a string or a type that can be converted to a string.
/// // let random_bool = next_random_bool!(seed_resource, "MyCategory");
/// ```
///
/// ## With the default category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`
/// // let random_bool = next_random_bool!(seed_resource);
/// ```
#[macro_export]
macro_rules! next_random_bool {
    // Variante für spezifische Kategorie
    ($seed_resource:expr, $category:expr) => {
        $seed_resource.next_bool(Some($category))
    };
    // Variante für Standard/Default-Wert (muss NACH der spezifischeren Variante stehen)
    ($seed_resource:expr) => {
        $seed_resource.next_bool(None)
    };
}

/// Generates the next random `bool` value based on a given probability for `true`, using the `SeedResource`.
///
/// The probability should be a `f32` value between `0.0` and `1.0`.
///
/// # Usage
///
/// ## With a specific category and probability:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`,
/// // `probability` is an f32, and `MyCategory` is a string.
/// // let random_bool = next_random_ratio!(seed_resource, probability, "MyCategory");
/// ```
///
/// ## With the default category and probability:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource` and `probability` is an f32.
/// // let random_bool = next_random_ratio!(seed_resource, probability);
/// ```
///
/// Consider renaming this macro to `next_random_bool_with_probability` for better clarity.
#[macro_export]
macro_rules! next_random_ratio {
    // Variante für spezifische Kategorie
    ($seed_resource:expr, $probability:expr, $category:expr) => {
        $seed_resource.next_bool_with_probability($probability, Some($category))
    };
    // Variante für Standard/Default-Wert (muss NACH der spezifischeren Variante stehen)
    ($seed_resource:expr, $probability:expr) => {
        $seed_resource.next_bool_with_probability($probability, None)
    };
}

/// Generates the next random `i64` integer within a specified range `[min, max)` (exclusive upper bound)
/// using the `SeedResource`.
///
/// # Usage
///
/// ## With a specific category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`,
/// // `min_val`, `max_val` are i64, and `MyCategory` is a string.
/// // let random_int = next_random_iter!(seed_resource, min_val, max_val, "MyCategory");
/// ```
///
/// ## With the default category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`,
/// // `min_val`, `max_val` are i64.
/// // let random_int = next_random_iter!(seed_resource, min_val, max_val);
/// ```
#[macro_export]
macro_rules! next_random_iter {
    // Exklusive Range: min..max
    ($seed_resource:expr, $min:expr, $max:expr, $category:expr) => {
        $seed_resource.next_i32_in_range($min, $max, Some($category)) // Assuming next_i32_in_range is the intended function
    };
    ($seed_resource:expr, $min:expr, $max:expr) => {
        $seed_resource.next_i32_in_range($min, $max, None) // Assuming next_i32_in_range is the intended function
    };
}

/// Generates the next random `f32` floating-point number within a specified range `[min, max)`
/// (exclusive upper bound, though typical for floats, precision might make it inclusive in some edge cases)
/// using the `SeedResource`.
///
/// # Usage
///
/// ## With a specific category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`,
/// // `min_val`, `max_val` are f32 or convertible to f32, and `MyCategory` is a string.
/// // let random_float = next_random_range!(seed_resource, min_val, max_val, "MyCategory");
/// ```
///
/// ## With the default category:
/// ```rust
/// // Assuming `seed_resource` is a mutable reference to a `SeedResource`,
/// // `min_val`, `max_val` are f32 or convertible to f32.
/// // let random_float = next_random_range!(seed_resource, min_val, max_val);
/// ```
#[macro_export]
macro_rules! next_random_range {
    // Exklusive Range mit Kategorie
    ($seed_resource:expr, $min:expr, $max:expr, $category:expr) => {
        $seed_resource.next_f32_in_range($min as f32, $max as f32, Some($category))
    };
    // Exklusive Range ohne Kategorie
    ($seed_resource:expr, $min:expr, $max:expr) => {
        $seed_resource.next_f32_in_range($min as f32, $max as f32, None)
    };
}
