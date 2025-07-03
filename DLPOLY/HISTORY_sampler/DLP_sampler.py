def extract_snapshots(input_file, output_file, interval):
    with open(input_file, 'r') as original_file:
        with open(output_file, 'w') as new_file:
            # Copying the first two lines of the original file to the new file
            for _ in range(2):
                line = original_file.readline()
                new_file.write(line)
                
            snapshot_count = 0
            while True:
                line = original_file.readline()
                if not line:
                    break
                if line.strip().startswith("timestep"):
                    snapshot_count += 1
                    if snapshot_count % interval == 1:
                        new_file.write(line)
                        next_line = original_file.readline()
                        while next_line and not next_line.strip().startswith("timestep"):
                            new_file.write(next_line)
                            next_line = original_file.readline()


input_file = "HISTORY"
output_file = "selected_snapshots_HISTORY"
interval = int(input("Enter the interval (n) for copying snapshots: "))

extract_snapshots(input_file, output_file, interval)
print("Selected snapshots with interval {} have been copied to {}.".format(interval, output_file))
