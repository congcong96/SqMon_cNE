# SqMon_cNE

The majority of codes are organized in **ne_analysis_main.m**  
The difference of this version from the 2018 paper:
  * NE events are identified as moments when more than one member neurons fire; 
  * Member neurona are labeled as positive or negative.

**sm_ne_NEmember_thresh_set** and **sm_ne_batch_save_NEtrain** are the mainfunctions that change the definition of NE events. The values that would have been obtained with 2018Paper method are still calculated and saved with suffix '_2018'
