# Load capture data and format 'captures' data frame for mrmr model

# Load .csv of all captures from the Mossy Pond CMR study period
captures <- read.csv("/Users/johnimperato/Documents/Capture_Data.csv")

# Remove "00:00:00 EDT" from Date column for compatibility with mrmr 
captures$Date <- gsub("00:00:00 EDT","", captures$Date)

# Convert 'Date' variable to date object class 
captures$Date <- as.Date(captures$Date, format = "%a %b %d  %Y")

# Change variable names for compatibility with mrmr 
names(captures)[names(captures) == "Date"] <- "survey_date"
names(captures)[names(captures) == "PIT_TagCode"] <- "pit_tag_id"

# Change data types for compatibility with mrmr
captures$survey_date <- as.character(captures$survey_date)
captures$pit_tag_id <- as.character(captures$pit_tag_id)

# Remove the single dead frog observation (instead, use it to create a 'removals' data frame)
captures <- captures[captures$frog_state != "dead", ]

# ATTEMPT TO FIGURE THE PROBLEM OUT
# Include only a subset of columns from captures. This seems to have worked in the past. 
### captures <- unique(captures[ , c(2,3,5,7)])
# when I create this subset, 5 observations disappear
# This subset didn't work. Lets try getting rid of CA lake id column
captures <- unique(captures[ , c(3,5)])
# another 2 observations disappeared
# I have not yet run this model

write.csv( captures , '/Users/johnimperato/Desktop/Mossy-Pond_companion/MP_captures.csv' )

