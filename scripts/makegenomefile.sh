bedtools makewindows -g genomefile.chr.txt -w 5000000 | awk '{print ":""-"}'|sed 's/:0-/:1-/' > genomefile.5M.txt
