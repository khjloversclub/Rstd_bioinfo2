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
spread(clinical_data, key = Treatment, value = Value)
patients$Name
patients[,c('Name', 'Sex')]
select(patients, -Name, Sex, Birth)
patients[,-2]
patients[,setdiff(colnames(patients), 'Name')]
select(patients, Name:Sex)
select(patients, starts_with('Grade'))
select(patients, Name, Sex:State, -Smokes)
patients <- read.delim('patient-data.txt')
patients <- as_tibble((patients))
patients
ggplot(patients, mapping=aes(x=Sex)) +geom_bar()
levels(patients$Sex)
patients$Sex <- str_trim(patients$Sex)
mutate(patients, Sex=str_trim(Sex))
mutate(patients, SexTrimmed=str_trim(Sex))
patients <- mutate(patients, Sex=factor(str_trim(Sex)))
select(patients, Height)
as.numeric(str_remove(patients$Height, 'cm'))
patients<- mutate(patients, 
                  Height=as.numeric(str_remove(patients$Height, 'cm')))
ggplot(patients, mapping=aes(x=Height)) +geom_histogram()
patients<- mutate(patients, 
                  Sex= as.factor(str_trim(Sex)),
                  Height=as.numeric(str_remove(Height, 'cm')))
patients<- read_tsv('patient-data-cleaned.txt')
mutate_at(patients, vars(Height, Weight), round, digits=1)
patients<- mutate_if(patients, is.numeric, round, digits=1)
patients


patients <- read_tsv("patient-data.txt")
patients %>% mutate(patients,Height=as.numeric(str_remove(Height, pattern='cm$')))
patients <- read.delim('patient-data.txt') %>%
  as_tibble %>%
  mutate(Sex= as_factor(str_trim(Sex))) %>%
  mutate(Height=as.numeric(str_remove(Height, pattern='cm$')))
patients
filter(patients, Sex == 'Male')
patients[patients$Sex == 'Male',]
filter(patients, Sex != 'Female')
filter(patients, State %in% c('Florida','Georgia','Illinois'))
filter(patients, str_detect(Name, '^B'))
filter(patients, str_starts(Name, 'B'))
filter(patients,!Died)
filter(patients, Sex=='Male'& Died)
patients[patients$Sex == 'Male' & patients$Died,]
filter(patients,Weight >90 | BMI >28, Sex=='Female')
arrange(patients, desc(Height))
arrange(patients,desc(Grade), Sex, Smokes)
candidates <- patients %>%
  filter(!Died) %>%
  select(ID,Name,Sex,Smokes,Height,Weight) %>%
  mutate(BMI = Weight / (Height / 100) ** 2) %>%
  mutate(Overweight = BMI > 25) %>%
  arrange(BMI)
candidates
