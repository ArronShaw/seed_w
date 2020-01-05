import os

ENTER_FOLDER = "Entering folder: "
#TEST_FOLDER = ["01_one_word", "02_two_words", "03_three_words"]
TEST_FOLDER = ["str_len", "str_num", "range", "word_num"]
TIME_OUTPUT = "time.out"
CSV_OUT = "results.csv"
# TIME_FORMAT = "-f '%C \nreal:%e \nuser:%U \nsys:%S \nMaxMem(Kb):%M \nAveMem(Kb):%t \nAveTolMem(Kb):%K' -a -o " + TIME_OUTPUT
TIME_FORMAT = "-f '%C \n%e \n%U \n%S \n%M \n%t \n%K' -a -o " + TIME_OUTPUT
# TIME_FORMAT = "-f '%C \n' -v -a -o " + TIME_OUTPUT

#str_len = ['1000', '2000', '4000', '8000', '16000', '32000', '64000']
str_len = ['1000']
#str_num = ['100', '200', '400', '800', '1600', '3200', '6400']
str_num = ['100']
# range_len = ['10', '20', '40', '80', '160', '320', '640']
range_len = ['640']
# word_num = ['1', '2', '3', '4', '6', '8']
word_num = ['1']

def print_pwd():
    print (ENTER_FOLDER + os.getcwd())

def excute(num):
    SEED_PARAMS = " --min_support 1" + " --word_length 10" \
    + " --mis_match_motif 2" + " --match_file save_motifs.xml"
    DATA_POS= " " + str(num) + "_pos.fas"
    DATA_NEG= " " + str(num) + "_neg.fas"
    # /usr/bin/time -v  on Linux
    # /time -l          on MacOS
    cmd_1 = "/usr/bin/time " + TIME_FORMAT 
    cmd_2 = " ../../seed_w " + SEED_PARAMS + DATA_POS + DATA_NEG + " > ./run.out"
    cmd = cmd_1 + cmd_2 
    #print(cmd)
    os.system(cmd)

def collect_data():
    
    with open(TIME_OUTPUT) as f:
        f.readline 
        for line in f:
            print(line, end='')
            
        f.close()
    print ("Results saved in:" + CSV_OUT)


# Enter folder recurcively and time the execution

#os.chdir("./examples_seed_w")
os.chdir("./src/algorithms/test_files")
print_pwd()
PARENT_FOLDER = os.getcwd()

for dir in TEST_FOLDER:
    os.chdir(PARENT_FOLDER + '/' + dir)
    print_pwd()
    # Clear TIME_OUTPUT File
    open(TIME_OUTPUT, "w").close()

    if dir == TEST_FOLDER[0]:
        len = str_len
    elif dir == TEST_FOLDER[1]:
        len = str_num
    elif dir == TEST_FOLDER[2]:
        len = range_len
    elif dir == TEST_FOLDER[3]:
        len = word_num
    else: 
        len = 0
    
    for i in len:
        excute(i)
        # break
    
    # break #test_folder
    #os.chdir(PARENT_FOLDER)


# Collect results
# os.chdir(PARENT_FOLDER)

# for dir in TEST_FOLDER:
#     os.chdir(dir)
#     print_pwd()
#     collect_data()
#     break