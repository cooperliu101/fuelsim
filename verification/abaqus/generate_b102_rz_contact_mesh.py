"""Create two independent QUAD8 rings with quadratic top and bottom contact edges."""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

root = Path(__file__).resolve().parent.parent / "meshes"
lower = [(1,0),(2,0),(2,1),(1,1),(1.5,0),(2,0.5),(1.5,1),(1,0.5)]
upper = [(1,1),(2,1),(2,2),(1,2),(1.5,1),(2,1.5),(1.5,2),(1,1.5)]
points = lower + upper
with netcdf_file(str(root / "b102_rz_contact.e"), "w") as f:
    f.api_version = np.float32(7.22); f.version = np.float32(7.22)
    f.floating_point_word_size = np.int32(8); f.file_size = np.int32(1)
    f.title = "CAX8T quadratic contact"
    for name, size in [("time_step",None),("len_name",33),("num_dim",2),("num_nodes",16),
                       ("num_elem",2),("num_el_blk",2),("num_node_sets",4),("num_side_sets",2),
                       ("num_el_in_blk1",1),("num_nod_per_el1",8),("num_el_in_blk2",1),
                       ("num_nod_per_el2",8),("num_side_ss1",1),("num_df_ss1",0),
                       ("num_side_ss2",1),("num_df_ss2",0)]:
        f.createDimension(name,size)
    def var(name,dims,values,kind="i"):
        v=f.createVariable(name,kind,dims);v[:]=values;return v
    def names(name,dim,values):
        out=np.zeros((len(values),33),dtype="S1")
        for i,value in enumerate(values): out[i,:len(value)]=np.frombuffer(value.encode(),dtype="S1")
        var(name,(dim,"len_name"),out,"c")
    var("coordx",("num_nodes",),[p[0] for p in points],"d");var("coordy",("num_nodes",),[p[1] for p in points],"d")
    var("connect1",("num_el_in_blk1","num_nod_per_el1"),[[1,2,3,4,5,6,7,8]]).elem_type="QUAD8"
    var("connect2",("num_el_in_blk2","num_nod_per_el2"),[[9,10,11,12,13,14,15,16]]).elem_type="QUAD8"
    var("eb_prop1",("num_el_blk",),[1,2]).name="ID";var("eb_status",("num_el_blk",),[1,1]);names("eb_names","num_el_blk",["lower","upper"])
    var("ss_prop1",("num_side_sets",),[1,2]).name="ID";var("ss_status",("num_side_sets",),[1,1]);names("ss_names","num_side_sets",["lower_contact","upper_contact"])
    var("elem_ss1",("num_side_ss1",),[1]);var("side_ss1",("num_side_ss1",),[3]);
    var("elem_ss2",("num_side_ss2",),[2]);var("side_ss2",("num_side_ss2",),[1]);
    var("ns_prop1",("num_node_sets",),[1,2,3,4]).name="ID";var("ns_status",("num_node_sets",),[1,1,1,1]);names("ns_names","num_node_sets",["lower_fixed","upper_fixed","lower_temp","upper_temp"])
    for n,nodes in [(1,[1,4]),(2,[11,12]),(3,[1,2,3,4,5,6,7,8]),(4,[9,10,11,12,13,14,15,16])]:
        dim="num_nod_ns%d"%n;f.createDimension(dim,len(nodes));var("node_ns%d"%n,(dim,),nodes)
