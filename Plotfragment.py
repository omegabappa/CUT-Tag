Use the below to find the fragment size distribution:
#!/usr/bin/env python3
import os
import pysam
import matplotlib.pyplot as plt

bam_dir = "sorted"
out_dir = "fragment_plots"
os.makedirs(out_dir, exist_ok=True)

for fn in sorted(os.listdir(bam_dir)):
    if not fn.endswith("_sorted.bam"):
        continue
    sample = fn.replace("_sorted.bam", "")
    bam_path = os.path.join(bam_dir, fn)
    sizes = []
    with pysam.AlignmentFile(bam_path, "rb") as bf:
        for r in bf.fetch():
            if r.is_proper_pair and r.template_length > 0:
                sizes.append(r.template_length)

    # plot
    plt.figure(figsize=(10, 6))
    plt.hist(sizes, bins=100)
    plt.title(f"Fragment Size Distribution — {sample}")
    plt.xlabel("Fragment Size (bp)")
    plt.ylabel("Count")
    plt.grid(True)
    plt.savefig(os.path.join(out_dir, f"{sample}_fragment_sizes.pdf"))
    plt.close()

##And run as:
##python3 plotfrag_python.py
