pub mod events;
pub mod plugin;
pub mod resource;

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

#[macro_export]
macro_rules! next_random_normalized {
    // Variante für spezifische Kategorie
    ($seed_resource_mut:expr, $category:expr) => {
        $seed_resource_mut.next_value_normalized($category)
    };
    // Variante für Standard/Default-Wert (muss NACH der spezifischeren Variante stehen)
    ($seed_resource_mut:expr) => {
        $seed_resource_mut.next_value_normalized()
    };
}

#[macro_export]
macro_rules! next_random_bool {
    // Variante für spezifische Kategorie
    ($seed_resource:expr, $category:expr) => {
        $seed_resource.next_bool($category)
    };
    // Variante für Standard/Default-Wert (muss NACH der spezifischeren Variante stehen)
    ($seed_resource:expr) => {
        $seed_resource.next_bool()
    };
}

#[macro_export]
macro_rules! next_random_ratio {
    // Variante für spezifische Kategorie
    ($seed_resource:expr, $category:expr) => {
        $seed_resource.next_bool_with_probability($category)
    };
    // Variante für Standard/Default-Wert (muss NACH der spezifischeren Variante stehen)
    ($seed_resource:expr) => {
        $seed_resource.next_bool_with_probability()
    };
}

#[macro_export]
macro_rules! next_random_iter {
    // Exklusive Range: min..max
    ($seed_resource:expr, $min:expr, $max:expr, $category:expr) => {
        $seed_resource.next_i64_in_range($min, $max, Some($category))
    };
    ($seed_resource:expr, $min:expr, $max:expr) => {
        $seed_resource.next_i64_in_range($min, $max, None)
    };
}

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
