patients <- read_tsv("patient-data-cleaned.txt")
View(patients)
library(tidyverse)
?geom_point
ggplot(data = patients) + geom_point(mapping= aes(x= BMI,y= Weight, colour = Smokes))
ggsave("newplot.png")

untidy_data <- read_tsv("tidyr-example.txt")
untidy_data
tidy_data <- gather(untidy_data, key="Treatment", value="Result",
                    treatmenta, treatmentb)
gather(untidy_data, Treatment, Result, treatmenta:treatmentb)
gather(untidy_data, Treatment, Result, -Name)
gather(untidy_data, Treatment, Result, 2:3)
clinical_data <- read_tsv("clinical-data.txt")
clinical_data <- gather(clinical_data, Treatment, Value, -Subject)
ggplot(data= clinical_data, mapping= aes(x=Treatment, y= Value))+ geom_boxplot()
clinical_replicate_data <- separate(clinical_data, Treatment, 
                                    into = c("Treatment", "Replicate"))
ggplot(data= clinical_replicate_data, 
       mapping= aes(x=Replicate, y= Value))+ geom_boxplot()+ facet_wrap(~Treatment)
