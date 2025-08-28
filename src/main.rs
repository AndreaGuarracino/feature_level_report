use clap::{Arg, Command};
use regex::Regex;
use std::{fs::File, io::{self, BufRead, BufReader}};
use flate2::read::GzDecoder;

fn count_aligned_bases(query_start: i64, query_end: i64, query_strand: char, target_start: i64, _target_end: i64, cigar: &str, feature_in_query_start: i64, feature_in_query_end: i64, feature_in_target_start: i64, feature_in_target_end: i64, max_indel_size: i64) -> (i64, i64, i64, i64, i64, i64, i64) {
    let mut aligned_bases = 0;
    let mut not_aligned_bases_in_query = 0;
    let mut not_aligned_bases_in_target = 0;
    let mut indels_in_query = 0;
    let mut indels_in_target =  0;
    let query_rev = query_strand == '-';

    // Initialize counters for the current position within the query and target sequences
    let mut query_pos = if query_rev { query_end } else { query_start };
    let mut target_pos = target_start;

    // Iterate over CIGAR operations
    let cigar_re = Regex::new(r"(\d+)([MIDNSHP=X])").unwrap();
    for cap in cigar_re.captures_iter(cigar) {
        let length = cap[1].parse::<i64>().unwrap();
        let op = &cap[2];
        match op {
            "M" | "=" | "X" => {
                // Handle match/mismatch, which affects both query and target
                let overlap_query = if query_rev {
                    std::cmp::max(0, std::cmp::min(query_pos, feature_in_query_end) - std::cmp::max(query_pos - length, feature_in_query_start))
                } else {
                    std::cmp::max(0, std::cmp::min(query_pos + length, feature_in_query_end) - std::cmp::max(query_pos, feature_in_query_start))
                };
                let overlap_target = std::cmp::max(0, std::cmp::min(target_pos + length, feature_in_target_end) - std::cmp::max(target_pos, feature_in_target_start));
                aligned_bases += std::cmp::min(overlap_query, overlap_target);

                if query_rev {
                    query_pos -= length;
                } else {
                    query_pos += length;
                }
                target_pos += length;
            },
            "D" => {
                // Handle deletion in the query (insertion in the target)
                let overlap_target = std::cmp::max(0, std::cmp::min(target_pos + length, feature_in_target_end) - std::cmp::max(target_pos, feature_in_target_start));
                if length <= max_indel_size {
                    indels_in_target += overlap_target;
                } else {
                    not_aligned_bases_in_target += overlap_target;
                }

                target_pos += length;
            },
            "I" => {
                // Handle insertion in the query (gap in the target)
                let overlap_query = if query_rev {
                    std::cmp::max(0, std::cmp::min(query_pos, feature_in_query_end) - std::cmp::max(query_pos - length, feature_in_query_start))
                } else {
                    std::cmp::max(0, std::cmp::min(query_pos + length, feature_in_query_end) - std::cmp::max(query_pos, feature_in_query_start))
                };
                if length <= max_indel_size {
                    indels_in_query += overlap_query;
                } else {
                    not_aligned_bases_in_query += overlap_query;
                }

                if query_rev {
                    query_pos -= length;
                } else {
                    query_pos += length;
                }
            },
            _ => {}
        }

        // Check if we have already passed the features in both query and target
        if ((query_rev && query_pos <= feature_in_query_start) || (!query_rev && query_pos >= feature_in_query_end)) && target_pos >= feature_in_target_end {
            break;
        }
    }
    (aligned_bases, not_aligned_bases_in_query, not_aligned_bases_in_target, indels_in_query, indels_in_target, (feature_in_query_end - feature_in_query_start) - aligned_bases - indels_in_query - not_aligned_bases_in_query, (feature_in_target_end - feature_in_target_start) - aligned_bases - indels_in_target - not_aligned_bases_in_target)
}

fn open_file(file_path: &str) -> Box<dyn BufRead> {
    if file_path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(File::open(file_path).expect("Failed to open file"))))
    } else {
        Box::new(BufReader::new(File::open(file_path).expect("Failed to open file")))
    }
}

fn main() -> io::Result<()> {
    let matches = Command::new("Alignment Feature Counter")
        .version("1.0")
        .author("Andrea Guarracino <aguarra1@uthsc.edu>")
        .about("Counts aligned bases for features in alignment data")
        .arg(Arg::new("input")
            .short('i')
            .long("input")
            .value_name("FILE")
            .help("Input file, can be gzipped")
            .num_args(1))
        .arg(Arg::new("max_indel_size")
            .short('m')
            .long("max-indel-size")
            .value_name("INT")
            .help("Maximum size of indels to consider in feature intervals")
            .num_args(1))
        .get_matches();

    let input_file = matches.get_one::<String>("input").map(|s| s.as_str()).unwrap_or("");
    let max_indel_size = matches.get_one::<String>("max_indel_size")
        .map(|s| s.parse::<i64>().expect("Invalid value for max indel size"))
        .unwrap_or(i64::MAX);

    println!("feature.name\tquery\tquery.feature.start\tquery.feature.end\tquery.strand\ttarget\ttarget.feature.start\ttarget.feature.end\taligned.bp\tnot.aligned.in.query.bp\tnot.aligned.in.target.bp\tindels.in.query.bp\tindels.in.target.bp\tignored.in.query.bp\tignored.in.target.bp");

    if !input_file.is_empty() {
        let file = open_file(input_file);
        for line in file.lines() {
            let line = line?;
            // Assuming `line` is a String obtained from iterating over lines of the file
            let parts: Vec<&str> = line.split('\t').collect();

            // Ensure there are enough parts to unpack
            // if parts.len() < 27 {
            //     eprintln!("ERROR: Line does not contain enough fields.");
            //     std::process::exit(1);
            // }

            let query_name = parts[0];
            //let query_len = parts[1].parse::<i64>().expect("Invalid query len");
            let query_start = parts[2].parse::<i64>().expect("Invalid query start");
            let query_end = parts[3].parse::<i64>().expect("Invalid query end");
            let query_strand = parts[4];
            let target_name = parts[5];
            //let target_len = parts[6].parse::<i64>().expect("Invalid target len");
            let target_start = parts[7].parse::<i64>().expect("Invalid target start");
            let target_end = parts[8].parse::<i64>().expect("Invalid target end");
            //_
            //_
            //_
            let cigar = parts[12].split("cg:Z:").last().unwrap_or_default();
            let query_name_2 = parts[13];
            let feature_in_query_start = parts[14].parse::<i64>().expect("Invalid feature in query start");
            let feature_in_query_end = parts[15].parse::<i64>().expect("Invalid feature in query end");
            let feature_in_query_name = parts[16];
            //_
            let feature_in_query_strand = parts[18];
            //let feature_in_query_class = parts[19];
            let target_name_2 = parts[20];
            let feature_in_target_start = parts[21].parse::<i64>().expect("Invalid feature in target start");
            let feature_in_target_end = parts[22].parse::<i64>().expect("Invalid feature in target end");
            let feature_in_target_name = parts[23];
            //_
            let feature_in_target_strand = parts[25];
            //let feature_in_target_class = parts[26];

            // Checking for matching names and strands
            if query_name != query_name_2 || target_name != target_name_2 || feature_in_query_name != feature_in_target_name {
                eprintln!("WARNING: query, target, and/or feature name do not match! Skip this line: {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", feature_in_query_name, query_name, feature_in_query_start, feature_in_query_end, query_strand, target_name, feature_in_target_start, feature_in_target_end);
                continue;
            }
            if feature_in_query_strand != feature_in_target_strand && query_strand == "+" {
                // If the features are on different strands, the query should be reversed in order to align them
                eprintln!("WARNING: the feature is on different strands in query and target, but query and target are in the same orientation! Skip this line: {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", feature_in_query_name, query_name, feature_in_query_start, feature_in_query_end, query_strand, target_name, feature_in_target_start, feature_in_target_end);
                continue;
            }

            let (aligned_bases, not_aligned_bases_in_query, not_aligned_bases_in_target, indels_in_query, indels_in_target, ignored_bases_in_query, ignored_bases_in_target) = count_aligned_bases(
                query_start, query_end, query_strand.chars().next().unwrap(), target_start, target_end, cigar, feature_in_query_start, feature_in_query_end, feature_in_target_start, feature_in_target_end, max_indel_size
            );

            println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", feature_in_query_name, query_name, feature_in_query_start, feature_in_query_end, query_strand, target_name, feature_in_target_start, feature_in_target_end, aligned_bases, not_aligned_bases_in_query, not_aligned_bases_in_target, indels_in_query, indels_in_target, ignored_bases_in_query, ignored_bases_in_target);
        }
    }

    Ok(())
}
