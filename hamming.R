# hamming.R
# Author : Mathieu G.
# Date : 3/12/2024
# Description : Function that calculates the hamming distance from 2 strings.
#               To do so, split strings into lists. Then compare the char at
#               each positions, add to a counter when there is a difference.
#               Return the counter.

# PRECONDITION : the strings must be of the same length.
# RETURNS : a number.

hamming <- function(str1, str2){
    
    count <- 0
    
    list_str1 <- unlist(strsplit(str1, ""))
    list_str2 <- unlist(strsplit(str2, ""))
    
    if (length(str1) != length(str2)) {
        return(NULL)  
    }
    
    else {
        for (i in 1:length(list_str1)) {
            if (list_str1[i] != list_str2[i]) {
                count <- count + 1
            }
        }
    }
    return(count)
}


