import os

VERSION = "1.0"

"""Output config"""
prefix = "gmw"
output_path = "./gmw_output"


"""threads config"""
max_threads = os.cpu_count()
#max_threads = 12

#overlap_size = 77

"""unfold position config"""
position_distance = 150
gc_discrepancy = 0.2
depth_discrepancy = 4

"""remove short reads offset"""
short_offset = 400

"""visualize config"""
edge_color="lavender"
node_shape="dot"
length_range_mix=50
length_range_max=500
width_range_mix=5
width_range_max=20
color_dict = {"unknown": "grey", "target": "pink", "contaminate": "skyblue"}

"""configs for kraken2"""
confidence=0.8
kraken_path = "kraken2"

"""blastn_config"""
blast_path = "blastn"
query_cover_offset=60
match_length=100