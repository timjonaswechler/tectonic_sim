// src/math/algorithms/convex_hull.rs

use crate::math::{
    error::{MathError, MathResult},
    utils::constants,
};
use bevy::math::Vec2;

/// Verschiedene Algorithmen zur Berechnung der konvexen Hülle einer Punktmenge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvexHullAlgorithm {
    GrahamScan,     // O(n log n)
    AndrewMonotone, // O(n log n), oft etwas schneller in der Praxis als Graham Scan
    JarvisMarch,    // O(nh), "Gift Wrapping", gut für wenige Hüllenpunkte (h)
    QuickHull,      // O(n log n) im Durchschnitt, O(n^2) im Worst-Case
}

impl Default for ConvexHullAlgorithm {
    fn default() -> Self {
        ConvexHullAlgorithm::AndrewMonotone // Guter Allrounder
    }
}

/// Berechnet die konvexe Hülle einer Menge von 2D-Punkten.
pub struct ConvexHullComputer {
    algorithm: ConvexHullAlgorithm,
    /// Ob kollineare Punkte auf dem Rand der Hülle mit eingeschlossen werden sollen.
    include_collinear_points: bool,
    /// Toleranz für Fließkommavergleiche, z.B. um Kollinearität zu bestimmen.
    tolerance: f32,
}

impl Default for ConvexHullComputer {
    fn default() -> Self {
        Self {
            algorithm: ConvexHullAlgorithm::default(),
            include_collinear_points: false,
            tolerance: constants::EPSILON * 10.0, // Etwas großzügiger für geometrische Vergleiche
        }
    }
}

impl ConvexHullComputer {
    pub fn new(algorithm: ConvexHullAlgorithm) -> Self {
        Self {
            algorithm,
            ..Default::default()
        }
    }

    pub fn include_collinear(mut self, include: bool) -> Self {
        self.include_collinear_points = include;
        self
    }

    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance.max(0.0);
        self
    }

    /// Berechnet die konvexe Hülle der gegebenen Punkte.
    /// Gibt eine Liste von Punkten zurück, die die Hülle in Gegen-Uhrzeigersinn-Orientierung bilden.
    /// Der erste und letzte Punkt sind nicht identisch in der zurückgegebenen Liste.
    /// Das Erstellen eines geschlossenen `Polygon`s daraus würde das Duplizieren des ersten Punktes erfordern.
    pub fn compute_hull_points(&self, input_points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        if input_points.len() < 3 {
            // Für weniger als 3 Punkte ist die "Hülle" die Menge der Punkte selbst,
            // oder man könnte einen Fehler werfen, da keine Fläche aufgespannt wird.
            // Hier geben wir die Punkte zurück, da sie eine Art "degenerierte" Hülle bilden.
            // Eine alternative wäre:
            // return Err(MathError::InsufficientPoints { expected: 3, actual: input_points.len() });
            // Für Konsistenz mit Polygon::closed, das 3 Punkte braucht, ist ein Fehler hier sinnvoll.
            if input_points.is_empty() {
                return Ok(Vec::new());
            } // Leere Eingabe, leere Ausgabe

            // Entferne Duplikate, bevor entschieden wird
            let mut unique_points = input_points.to_vec();
            unique_points.sort_by(|a, b| {
                a.x.partial_cmp(&b.x)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
            });
            unique_points.dedup_by(|a, b| a.distance_squared(*b) < self.tolerance * self.tolerance);

            if unique_points.len() < 3 {
                return Err(MathError::InsufficientPoints {
                    expected: 3,
                    actual: unique_points.len(),
                });
            }
            // Wenn nach Dedup genug Punkte da sind, trotzdem den Algorithmus laufen lassen.
            // Ansonsten ist die Hülle einfach die Menge der unique_points, wenn < 3.
            // Der Algorithmus selbst sollte das aber korrekt handhaben oder einen Fehler werfen.
        }

        // Entferne Duplikate vor der Verarbeitung, um Algorithmen zu vereinfachen
        let mut points = input_points.to_vec();
        points.sort_by(|a, b| {
            // Stabile Sortierung für konsistente Ergebnisse
            a.x.partial_cmp(&b.x)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
        });
        points.dedup_by(|a, b| a.distance_squared(*b) < self.tolerance * self.tolerance);

        if points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: points.len(),
            });
        }

        let hull_points = match self.algorithm {
            ConvexHullAlgorithm::GrahamScan => self.graham_scan(&points)?,
            ConvexHullAlgorithm::AndrewMonotone => self.andrew_monotone(&points)?,
            ConvexHullAlgorithm::JarvisMarch => self.jarvis_march(&points)?,
            ConvexHullAlgorithm::QuickHull => self.quick_hull(&points)?,
        };

        if hull_points.len() < 3 {
            // Sollte nicht passieren, wenn Eingabe > 3 war und Algorithmus korrekt ist,
            // aber als Sicherheitsnetz.
            return Err(MathError::GeometricFailure {
                operation: "Convex hull computation resulted in fewer than 3 unique points."
                    .to_string(),
            });
        }
        Ok(hull_points)
    }

    // --- Private Hilfsfunktion: Kreuzprodukt (orientierter Flächeninhalt) ---
    // Gibt > 0 zurück, wenn c links von der Geraden ab liegt (CCW-Drehung)
    // Gibt < 0 zurück, wenn c rechts von der Geraden ab liegt (CW-Drehung)
    // Gibt ~0 zurück, wenn a, b, c kollinear sind.
    #[inline]
    fn cross_product_orientation(p: Vec2, q: Vec2, r: Vec2) -> f32 {
        (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
    }

    // --- Graham Scan ---
    fn graham_scan(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        let mut local_points = points.to_vec(); // Kopie zum Sortieren
        let n = local_points.len();

        // 1. Finde den Punkt p0 mit der kleinsten y-Koordinate.
        //    Bei Gleichheit den mit der kleinsten x-Koordinate.
        let mut min_idx = 0;
        for i in 1..n {
            if local_points[i].y < local_points[min_idx].y
                || (local_points[i].y == local_points[min_idx].y
                    && local_points[i].x < local_points[min_idx].x)
            {
                min_idx = i;
            }
        }
        local_points.swap(0, min_idx);
        let p0 = local_points[0];

        // 2. Sortiere die restlichen Punkte nach Polarwinkel relativ zu p0.
        //    Bei gleichem Winkel, behalte den fernsten Punkt (oder den nächsten, je nach include_collinear).
        local_points[1..].sort_unstable_by(|pa, pb| {
            let orientation = Self::cross_product_orientation(p0, *pa, *pb);
            if orientation.abs() < self.tolerance {
                // Kollinear
                // Behalte den fernsten Punkt, wenn include_collinear false ist, oder den nächsten, wenn true.
                // Standardmäßig (false): fernster Punkt bevorzugt, um Kollinearität auf der Hülle zu reduzieren.
                let dist_sq_a = p0.distance_squared(*pa);
                let dist_sq_b = p0.distance_squared(*pb);
                if self.include_collinear_points {
                    dist_sq_a
                        .partial_cmp(&dist_sq_b)
                        .unwrap_or(std::cmp::Ordering::Equal)
                } else {
                    dist_sq_b
                        .partial_cmp(&dist_sq_a)
                        .unwrap_or(std::cmp::Ordering::Equal) // umgekehrter Vergleich für fernsten
                }
            } else if orientation > 0.0 {
                // pa ist links von pb (CCW)
                std::cmp::Ordering::Less
            } else {
                // pa ist rechts von pb (CW)
                std::cmp::Ordering::Greater
            }
        });

        // Optional: Entferne Punkte, die nach der Sortierung bei gleichem Winkel näher an p0 sind,
        // falls include_collinear_points false ist. Die Sortierung sollte das schon erledigt haben.

        // 3. Baue die Hülle auf
        let mut hull: Vec<Vec2> = Vec::with_capacity(n);
        hull.push(p0);
        if n > 1 {
            hull.push(local_points[1]);
        }

        for i in 2..n {
            let mut top = hull.pop().unwrap();
            // Solange der letzte Turn keine Links-Drehung ist (oder eine tolerierte Rechts-Drehung/Kollinearität)
            while !hull.is_empty() {
                let orientation = Self::cross_product_orientation(
                    hull.last().unwrap().clone(),
                    top,
                    local_points[i],
                );
                if orientation > self.tolerance {
                    // Strikt links
                    break;
                }
                if !self.include_collinear_points && orientation.abs() < self.tolerance {
                    // Kollinear und nicht erwünscht
                    // Ersetze `top` mit dem weiter entfernten Punkt (current_points[i]), wenn es weiter von hull.last() ist
                    if hull.last().unwrap().distance_squared(local_points[i])
                        > hull.last().unwrap().distance_squared(top)
                    {
                        // Diese Logik ist etwas anders als Standard Graham, der Kollinearität meist ignoriert oder
                        // einfach den letzten kollinearen Punkt behält.
                        // Hier versuchen wir, den äußersten kollinearen Punkt zu behalten.
                        // Standard Graham: Wenn kollinear, pop. Wenn include_collinear, nicht poppen.
                        // Da wir oben schon nach Distanz sortiert haben (fernster bei !include_collinear),
                        // sollte das hier einfacher sein.
                        top = hull.pop().unwrap(); // Standard Graham: poppen bei <=0
                        continue;
                    } else {
                        break; // Behalte `top`, da `local_points[i]` nicht weiter entfernt ist
                    }
                }
                if orientation < -self.tolerance {
                    // Eindeutig rechts
                    top = hull.pop().unwrap();
                    continue;
                }
                // Fall: Kollinear und include_collinear_points = true
                break;
            }
            hull.push(top);
            hull.push(local_points[i]);
        }
        Ok(hull)
    }

    // --- Andrew's Monotone Chain ---
    fn andrew_monotone(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        // Punkte sind bereits sortiert und dedupliziert durch compute_hull_points
        let n = points.len();

        let mut upper_hull: Vec<Vec2> = Vec::with_capacity(n);
        let mut lower_hull: Vec<Vec2> = Vec::with_capacity(n);

        for &p_i in points {
            // Untere Hülle
            while lower_hull.len() >= 2 {
                let orientation = Self::cross_product_orientation(
                    lower_hull[lower_hull.len() - 2],
                    lower_hull[lower_hull.len() - 1],
                    p_i,
                );
                if orientation < -self.tolerance
                    || (!self.include_collinear_points && orientation.abs() < self.tolerance)
                {
                    lower_hull.pop();
                } else {
                    break;
                }
            }
            lower_hull.push(p_i);

            // Obere Hülle
            while upper_hull.len() >= 2 {
                let orientation = Self::cross_product_orientation(
                    upper_hull[upper_hull.len() - 2],
                    upper_hull[upper_hull.len() - 1],
                    p_i,
                );
                // Für obere Hülle umgekehrte Logik: rechts-Turn oder Kollinearität (wenn nicht erlaubt) führt zum Pop
                if orientation > self.tolerance
                    || (!self.include_collinear_points && orientation.abs() < self.tolerance)
                {
                    upper_hull.pop();
                } else {
                    break;
                }
            }
            upper_hull.push(p_i);
        }

        // Kombiniere: untere Hülle (ohne letzten Punkt) + umgekehrte obere Hülle (ohne ersten und letzten Punkt)
        // Andrew's Standard: lower_hull + upper_hull.reverse()[1..n-1]
        let mut hull = lower_hull;
        // Füge die obere Hülle (außer dem ersten und letzten Punkt, da sie in lower_hull sind) in umgekehrter Reihenfolge hinzu.
        // Store the first and last points before extending
        let first_point = hull.first().copied();
        let last_point = hull.last().copied();

        hull.extend(
            upper_hull
                .into_iter()
                .rev()
                .skip(1)
                .take_while(|p| {
                    last_point.map_or(true, |lp| {
                        lp.distance_squared(*p) > self.tolerance * self.tolerance
                    })
                })
                .skip_while(|p| {
                    first_point.map_or(false, |fp| {
                        fp.distance_squared(*p) < self.tolerance * self.tolerance
                    })
                }),
        );
        // Die obige extend Logik ist kompliziert. Einfacher:
        // hull.extend_from_slice(&upper_hull[1..upper_hull.len()-1].iter().rev().cloned().collect::<Vec<_>>());
        // Noch einfacher (Standardweg):
        // let mut hull = lower_hull;
        // hull.pop(); // Entferne letzten Punkt von lower (ist Startpunkt von upper)
        // upper_hull.pop(); // Entferne letzten Punkt von upper (ist Endpunkt von lower, wenn umgedreht)
        // hull.extend(upper_hull.into_iter().rev());

        // Korrekter Weg für Andrew's:
        // 1. Build lower hull
        // 2. Build upper hull (iterating in reverse order of points)
        // 3. Concatenate lower (without its last point) and upper (without its last point)

        // Neuaufbau Andrew's Monotone Chain nach Standard-Algorithmus:
        let mut hull_final: Vec<Vec2> = Vec::new();
        // Untere Hülle
        for &p in points {
            while hull_final.len() >= 2 {
                let cross = Self::cross_product_orientation(
                    hull_final[hull_final.len() - 2],
                    hull_final[hull_final.len() - 1],
                    p,
                );
                if cross < -self.tolerance
                    || (!self.include_collinear_points && cross.abs() < self.tolerance)
                {
                    hull_final.pop();
                } else {
                    break;
                }
            }
            hull_final.push(p);
        }
        // Obere Hülle
        // Speichere die Größe der unteren Hülle, bevor wir die obere hinzufügen (außer dem letzten Punkt der unteren Hülle)
        let lower_hull_size = hull_final.len();
        // Iteriere rückwärts für die obere Hülle, beginnend beim vorletzten Punkt der sortierten Liste
        for i in (0..n - 1).rev() {
            // n-1 ist der letzte Punkt, n-2 ist der vorletzte
            let p = points[i];
            while hull_final.len() > lower_hull_size {
                // Stelle sicher, dass wir nicht die untere Hülle weiter reduzieren
                let cross = Self::cross_product_orientation(
                    hull_final[hull_final.len() - 2],
                    hull_final[hull_final.len() - 1],
                    p,
                );
                if cross < -self.tolerance
                    || (!self.include_collinear_points && cross.abs() < self.tolerance)
                {
                    hull_final.pop();
                } else {
                    break;
                }
            }
            hull_final.push(p);
        }
        // Der letzte Punkt der oberen Hülle (der ursprüngliche Startpunkt) ist jetzt doppelt, entferne ihn.
        if hull_final.len() > 1 && hull_final.first() == hull_final.last() {
            // Nur wenn wirklich dupliziert
            hull_final.pop();
        }

        Ok(hull_final)
    }

    // --- Jarvis March (Gift Wrapping) ---
    fn jarvis_march(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        let n = points.len();
        // Finde den Punkt mit der kleinsten y-Koordinate (linkester bei Gleichheit)
        let mut start_idx = 0;
        for i in 1..n {
            if points[i].y < points[start_idx].y
                || (points[i].y == points[start_idx].y && points[i].x < points[start_idx].x)
            {
                start_idx = i;
            }
        }

        let mut hull_indices = Vec::new();
        let mut current_hull_idx = start_idx;
        loop {
            hull_indices.push(current_hull_idx);
            let mut next_candidate_idx = (current_hull_idx + 1) % n; // Wähle irgendeinen anderen Punkt als ersten Kandidaten

            for i in 0..n {
                if i == current_hull_idx {
                    continue;
                }

                let orientation = Self::cross_product_orientation(
                    points[current_hull_idx],
                    points[i], // Teste diesen Punkt
                    points[next_candidate_idx],
                );

                if orientation < -self.tolerance {
                    // points[i] ist weiter rechts (CW) als next_candidate_idx
                    next_candidate_idx = i;
                } else if orientation.abs() < self.tolerance {
                    // Kollinear
                    // Wähle den Punkt, der weiter entfernt ist vom current_hull_idx
                    let dist_sq_i = points[current_hull_idx].distance_squared(points[i]);
                    let dist_sq_cand =
                        points[current_hull_idx].distance_squared(points[next_candidate_idx]);
                    if self.include_collinear_points {
                        if dist_sq_i > dist_sq_cand {
                            next_candidate_idx = i;
                        }
                    } else {
                        // Wenn nicht include_collinear, und es ist ein "echter" Turn (nicht zurück)
                        if dist_sq_i > dist_sq_cand {
                            next_candidate_idx = i;
                        }
                    }
                }
            }
            current_hull_idx = next_candidate_idx;
            if current_hull_idx == start_idx {
                break;
            } // Zurück am Startpunkt
            if hull_indices.len() > n {
                // Sicherheitsabbruch
                return Err(MathError::GeometricFailure {
                    operation: "Jarvis March did not terminate".to_string(),
                });
            }
        }
        Ok(hull_indices.into_iter().map(|idx| points[idx]).collect())
    }

    // --- QuickHull ---
    fn quick_hull(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        // Punkte sind bereits sortiert und dedupliziert
        let n = points.len();

        // Finde die Punkte mit minimaler und maximaler x-Koordinate
        // Da `points` nach x sortiert ist, sind das points[0] und points[n-1]
        let p_min_x = points[0];
        let p_max_x = points[n - 1];

        let mut hull = Vec::new();
        hull.push(p_min_x); // Start der Hülle

        // Teile die Punkte in zwei Mengen S1 und S2
        // S1: Punkte rechts von der Geraden (p_min_x -> p_max_x)
        // S2: Punkte rechts von der Geraden (p_max_x -> p_min_x) (also links von p_min_x -> p_max_x)
        let mut s1 = Vec::new();
        let mut s2 = Vec::new();

        for i in 1..(n - 1) {
            // Überspringe p_min_x und p_max_x
            let p = points[i];
            let orientation = Self::cross_product_orientation(p_min_x, p_max_x, p);
            if orientation > self.tolerance {
                // Links von (p_min_x -> p_max_x) -> Kandidat für obere Hülle
                s1.push(p);
            } else if orientation < -self.tolerance {
                // Rechts von (p_min_x -> p_max_x) -> Kandidat für untere Hülle
                s2.push(p);
            } else if self.include_collinear_points {
                // Behandle kollineare Punkte hier, wenn sie Teil der Hülle sein sollen.
                // Füge sie zu s1 oder s2 hinzu, je nachdem, ob sie zwischen p_min_x und p_max_x liegen.
                if p.x > p_min_x.x && p.x < p_max_x.x { // Nur wenn strikt dazwischen
                    // Hier ist es schwierig zu entscheiden, zu welcher Seite sie gehören, wenn sie exakt auf der Linie liegen.
                    // Man könnte sie ignorieren oder einer Seite zuordnen.
                    // Für include_collinear ist es besser, sie in der rekursiven Suche zu berücksichtigen.
                }
            }
        }

        self.quick_hull_recursive(&mut hull, &s1, p_min_x, p_max_x);
        hull.push(p_max_x); // Mittlerer Punkt der Hülle
        self.quick_hull_recursive(&mut hull, &s2, p_max_x, p_min_x); // Für die andere Seite

        // Der letzte Punkt könnte der erste sein, wenn die Hülle geschlossen ist. Polygon::closed kümmert sich darum.
        // Hier entfernen wir das Duplikat, falls es durch die Rekursion entstanden ist.
        if hull.len() > 1 && hull.first() == hull.last() {
            hull.pop();
        }

        Ok(hull)
    }

    fn quick_hull_recursive(
        &self,
        hull_segment: &mut Vec<Vec2>,
        points_set: &[Vec2],
        p_a: Vec2,
        p_b: Vec2,
    ) {
        if points_set.is_empty() {
            return;
        }

        // Finde den Punkt p_c im Set, der den größten Abstand zur Geraden (p_a, p_b) hat.
        // (oder den mit dem größten vorzeichenbehafteten Kreuzprodukt, um auf der richtigen Seite zu bleiben)
        let mut max_dist = -1.0; // Negativ, da wir nur Punkte auf einer Seite betrachten
        let mut p_c_idx = usize::MAX;

        for (idx, &p_candidate) in points_set.iter().enumerate() {
            // Wir erwarten, dass alle Punkte in points_set auf derselben Seite von (p_a, p_b) liegen.
            // Das Kreuzprodukt gibt uns einen Hinweis auf den Abstand und die Seite.
            let current_orientation_dist = Self::cross_product_orientation(p_a, p_b, p_candidate);
            // Quickhull sucht Punkte auf der "linken" Seite der gerichteten Kante p_a->p_b
            // D.h. cross_product_orientation > 0
            if current_orientation_dist > max_dist + self.tolerance {
                // +tolerance um den größten positiven zu finden
                max_dist = current_orientation_dist;
                p_c_idx = idx;
            } else if self.include_collinear_points
                && current_orientation_dist.abs() < self.tolerance
            {
                // Wenn kollinear und erlaubt, wähle den entferntesten von p_a
                if p_c_idx == usize::MAX
                    || p_a.distance_squared(p_candidate) > p_a.distance_squared(points_set[p_c_idx])
                {
                    max_dist = current_orientation_dist; // Bleibt nahe 0
                    p_c_idx = idx;
                }
            }
        }

        if p_c_idx == usize::MAX {
            // Kein Punkt links von der Linie gefunden
            return;
        }
        let p_c = points_set[p_c_idx];

        // Punkte rechts von (p_a, p_c) -> s_ac
        let mut s_ac = Vec::new();
        for &p_s in points_set {
            if Self::cross_product_orientation(p_a, p_c, p_s) > self.tolerance {
                s_ac.push(p_s);
            } else if self.include_collinear_points
                && Self::cross_product_orientation(p_a, p_c, p_s).abs() < self.tolerance
            {
                // Behandle Kollinearität für (p_a, p_c)
                if p_s.x > p_a.x.min(p_c.x)
                    && p_s.x < p_a.x.max(p_c.x)
                    && p_s.y > p_a.y.min(p_c.y)
                    && p_s.y < p_a.y.max(p_c.y)
                {
                    // Nur wenn dazwischen
                    s_ac.push(p_s);
                }
            }
        }
        self.quick_hull_recursive(hull_segment, &s_ac, p_a, p_c);

        hull_segment.push(p_c); // p_c ist Teil der Hülle

        // Punkte rechts von (p_c, p_b) -> s_cb
        let mut s_cb = Vec::new();
        for &p_s in points_set {
            if Self::cross_product_orientation(p_c, p_b, p_s) > self.tolerance {
                s_cb.push(p_s);
            } else if self.include_collinear_points
                && Self::cross_product_orientation(p_c, p_b, p_s).abs() < self.tolerance
            {
                // Behandle Kollinearität für (p_c, p_b)
                if p_s.x > p_c.x.min(p_b.x)
                    && p_s.x < p_c.x.max(p_b.x)
                    && p_s.y > p_c.y.min(p_b.y)
                    && p_s.y < p_c.y.max(p_b.y)
                {
                    s_cb.push(p_s);
                }
            }
        }
        self.quick_hull_recursive(hull_segment, &s_cb, p_c, p_b);
    }
}
