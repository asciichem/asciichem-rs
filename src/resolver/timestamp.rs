//! UTC timestamps for provenance and cache entries — the same
//! second-precision `%Y-%m-%dT%H:%M:%SZ` form the Ruby resolver
//! stamps (`Time.now.utc.strftime`) and parses (`Time.parse`, loose
//! on garbage). Plain civil-date arithmetic, no clock dependency:
//! the dates live in a small, well-tested algorithm (Howard
//! Hinnant's civil-from-days / days-from-civil).

use std::time::{SystemTime, UNIX_EPOCH};

/// The current UTC instant, second precision, wire form.
pub(crate) fn now_utc() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|since| since.as_secs() as i64)
        .unwrap_or_default();
    format_utc(secs)
}

/// Formats a unix-epoch second count as `YYYY-MM-DDTHH:MM:SSZ`.
pub(crate) fn format_utc(secs: i64) -> String {
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        rem / 3_600,
        (rem % 3_600) / 60,
        rem % 60
    )
}

/// Parses the wire timestamp form back to epoch seconds. Anything
/// that is not exactly the wire shape parses as nothing — the cache
/// treats that as stale.
pub(crate) fn parse_utc(text: &str) -> Option<i64> {
    let bytes = text.as_bytes();
    if text.len() != 20
        || bytes[4] != b'-'
        || bytes[7] != b'-'
        || bytes[10] != b'T'
        || bytes[13] != b':'
        || bytes[16] != b':'
        || bytes[19] != b'Z'
    {
        return None;
    }
    let (year, month, day) = (num(&text[0..4])?, num(&text[5..7])?, num(&text[8..10])?);
    let (hour, minute, second) = (
        num(&text[11..13])?,
        num(&text[14..16])?,
        num(&text[17..19])?,
    );
    if !(1..=12).contains(&month)
        || !(1..=31).contains(&day)
        || hour > 23
        || minute > 59
        || second > 59
    {
        return None;
    }
    Some(
        days_from_civil(year, month as u32, day as u32) * 86_400
            + hour * 3_600
            + minute * 60
            + second,
    )
}

fn num(text: &str) -> Option<i64> {
    text.parse::<i64>().ok()
}

/// Days since 1970-01-01 for a civil date (Hinnant's algorithm).
fn days_from_civil(year: i64, month: u32, day: u32) -> i64 {
    let y = if month <= 2 { year - 1 } else { year };
    let m = if month <= 2 { month + 9 } else { month - 3 } as i64;
    let era = y.div_euclid(400);
    let yoe = y - era * 400; // [0, 399]
    let doy = (153 * m + 2) / 5 + day as i64 - 1; // [0, 365]
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy; // [0, 146096]
    era * 146_097 + doe - 719_468
}

/// The civil date for a days-since-epoch count (Hinnant's algorithm).
fn civil_from_days(days: i64) -> (i64, u32, u32) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097; // [0, 146096]
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365; // [0, 399]
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32; // [1, 12]
    (if m <= 2 { y + 1 } else { y }, m, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_known_instants() {
        assert_eq!(format_utc(0), "1970-01-01T00:00:00Z");
        assert_eq!(format_utc(1_000_000_000), "2001-09-09T01:46:40Z");
        assert_eq!(format_utc(2_000_000_000), "2033-05-18T03:33:20Z");
    }

    #[test]
    fn round_trips_the_wire_form() {
        for secs in [0, 86_399, 86_400, 951_782_400, 1_762_000_000, 4_102_444_800] {
            assert_eq!(parse_utc(&format_utc(secs)), Some(secs), "{secs}");
        }
    }

    #[test]
    fn rejects_garbage_as_stale() {
        assert_eq!(parse_utc(""), None);
        assert_eq!(parse_utc("not a timestamp"), None);
        assert_eq!(parse_utc("2026-11-12 00:00:00"), None);
        assert_eq!(parse_utc("2026-13-40T99:99:99Z"), None);
    }
}
