library(ggplot2)
library(stringr)
library(cowplot)


ref <- readDNAStringSet("rawdata/G_wt.fasta")
refaa <- translate(ref)



# Beta strandi	13 – 17	Combined sources		5
# Beta strandi	26 – 28	Combined sources		3
# Beta strandi	30 – 33	Combined sources		4
# Beta strandi	36 – 49	Combined sources		14
# Beta strandi	51 – 61	Combined sources		11
# Beta strandi	69 – 83	Combined sources		15
# Beta strandi	88 – 98	Combined sources		11
# Beta strandi	108 – 110	Combined sources		3
# Beta strandi	115 – 117	Combined sources		3
# Beta strandi	120 – 130	Combined sources		11
# Beta strandi	134 – 137	Combined sources		4
# Beta strandi	139 – 150	Combined sources		12
# Beta strandi	152 – 163	Combined sources		12