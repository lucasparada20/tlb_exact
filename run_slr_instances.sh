
# Model can be 1,2,3,5,6,7
build/exec_tlb mcfp_solver=3 instance_format=trips \
    instance_file=some/path/to/montreal_n801_e15_0.txt \
    scenarios_file_name=some/path/to/montreal_e400.txt \
    model=1 \
    re_file=re.txt \
    output_file=targets.txt

build/exec_tlb mcfp_solver=3 instance_format=trips \
    instance_file=some/path/to/boston_n424_e15_0.txt \
    scenarios_file_name=some/path/to/boston_e400.txt \
    model=3 \
    re_file=re.txt \
    output_file=targets.txt

build/exec_tlb mcfp_solver=3 instance_format=trips \
    instance_file=some/path/to/newyork_n2176_e15_0.txt \
    scenarios_file_name=some/path/to/newyork_e400.txt \
    model=3 \
    re_file=re.txt \
    output_file=targets.txt
