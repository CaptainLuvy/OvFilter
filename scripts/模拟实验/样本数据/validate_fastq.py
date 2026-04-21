# -*- coding: utf-8 -*-
import os

# Use relative paths to avoid encoding issues with Chinese characters in the script
FILES_TO_CHECK = [
    "sample_profile_hifi.fastq",
    "sample_profile_ont.fastq"
]

def validate_and_fix_fastq(filename):
    file_path = os.path.abspath(filename)
    print(f"\nChecking: {filename}")
    
    if not os.path.exists(file_path):
        print("  Error: File not found.")
        return

    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"  Error reading file: {e}")
        return

    total_lines = len(lines)
    print(f"  Total lines: {total_lines}")

    if total_lines == 0:
        print("  Error: File is empty.")
        return

    # 1. Check for incomplete records at the end
    remainder = total_lines % 4
    if remainder != 0:
        print(f"  Warning: File has incomplete record at the end (remainder {remainder}).")
        print("  Fixing: Removing incomplete lines...")
        lines = lines[:-remainder]
        total_lines = len(lines)
        
        # Write back fixed content
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.writelines(lines)
            print(f"  Fixed. New line count: {total_lines}")
        except Exception as e:
            print(f"  Error writing file: {e}")
            return
    else:
        print("  Structure: Line count is valid (multiple of 4).")

    # 2. Check format of first and last record
    if total_lines >= 4:
        # First record check
        if not lines[0].startswith('@'):
            print("  Error: Line 1 does not start with '@'")
        elif not lines[2].startswith('+'):
            print("  Error: Line 3 does not start with '+'")
        else:
            # Check lengths for first record
            seq_len = len(lines[1].strip())
            qual_len = len(lines[3].strip())
            if seq_len != qual_len:
                print(f"  Error: First record sequence length ({seq_len}) != quality length ({qual_len})")
            else:
                print("  Content: First record format looks valid.")
        
        # Last record check
        last_idx = total_lines - 4
        if not lines[last_idx].startswith('@'):
            print(f"  Error: Line {last_idx+1} does not start with '@'")
        if not lines[last_idx+2].startswith('+'):
            print(f"  Error: Line {last_idx+3} does not start with '+'")

    print("  Status: Ready for use.")

if __name__ == "__main__":
    for f in FILES_TO_CHECK:
        validate_and_fix_fastq(f)
