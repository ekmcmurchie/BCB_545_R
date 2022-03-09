library(tidyverse)
dvst <- read_csv("https://raw.githubusercontent.com/vsbuffalo/bds-files/master/chapter-08-r/Dataset_S1.txt")
dvst
view(dvst)
read_csv()

summary(select(dvst,`total SNPs`)) # Here we also used `select` function. We'll talk about it soon.
filter(dvst,`total SNPs` >= 85)
filter(dvst, Pi > 16, `%GC` > 80)
(new_df <- filter(dvst, Pi > 16, `%GC` > 80))
mutate(dvst, cent = start >= 25800000 & end <= 29700000)
dvst <- mutate(dvst, 
               diversity = Pi / (10*1000), 
               cent = start >= 25800000 & end <= 29700000
)
view(dvst)
summary(select(dvst, `%GC`)) # lets us see highest and lowest %GC
filter(dvst, SNPs == 0)
mutate(dvst, no_SNPs = SNPs == 0)
view(dvst)

#ggplot2 tutorial
dvst <- read_csv("https://raw.githubusercontent.com/vsbuffalo/bds-files/master/chapter-08-r/Dataset_S1.txt") %>% 
  mutate(diversity = Pi / (10*1000), cent = (start >= 25800000 & end <= 29700000)) %>% 
  rename(percent.GC = `%GC`, total.SNPs = `total SNPs`, total.Bases = `total Bases`, reference.Bases = `reference Bases`)

dvst <- mutate(dvst, position = (end + start) / 2) #add position column
ggplot(data = dvst) + 
  geom_point(mapping = aes(x = position, y = diversity))
#plotting diversity vs. position
ggplot(data = dvst) +
  geom_point(mapping = aes(x = position, y = diversity, color = cent))
#color based on cetromere or not. hard to sequence centromere
ggplot(data = dvst) + 
  geom_point(mapping = aes(x = position, y = diversity), color = "blue")
#same thing but blue
ggplot(data = dvst, mapping = (aes(x = position, y = diversity))) + 
  geom_point()
#calling mapping first - good if using same mapping for different geom


ggplot(data = dvst) +
  geom_point(mapping = aes(x = position, y = diversity, color = start, size = start))
# can't set continuous variable for shape but can set size and color
?geom_point

ggplot(data = dvst) +
  geom_point(mapping = aes(x = position, y = diversity, color = cent), shape = 3, stroke = 3)

ggplot(data = dvst) +
  geom_point(mapping = aes(x = position, y = diversity), alpha = 0.01)
# alpha value adds transparency, to avoid overplotting

ggplot(data = dvst) +
  geom_density(mapping = aes(x = diversity), fill = "blue")
# density of diversity values

ggplot(data = dvst) +
  geom_density(mapping = aes(x = diversity, fill = cent), alpha = 0.4)
#fill color based on if in centromere or not

ggplot(data = dvst, mapping = aes(x = depth, y = total.SNPs)) +
  geom_point(alpha = 0.1) +
  geom_smooth()
# two geoms, both point and smoothing line

ggplot(data = dvst, mapping = aes(x = percent.GC, y = depth)) +
  geom_point(alpha = 0.1) +
  geom_smooth()

ggplot(dvst) +
  geom_point(aes(x = position, y = diversity)) +
  xlab("chromosome position (basepairs)") +
  ylab("nucleotide diversity")
# adding labels

ggplot(data = dvst) + 
  geom_bar(mapping = aes(x = cent))
# bar plot - centromere or not (ugly)

dvst <- dvst %>% 
  mutate(GC.binned = cut(percent.GC, 5));
select(dvst, GC.binned)
#here we cut percent.GC into 5 different bins, low to high, in a new column
#cuts continuous variable into categories
ggplot(data = dvst) + geom_bar(mapping = aes(x = GC.binned))

ggplot(data = dvst) +
  geom_density(mapping = aes(x = depth, linetype = GC.binned), alpha = 0.5)
# not a useful graph

ggplot(data = dvst) + geom_bar(mapping = aes(x = percent.GC))
#automatic histogram - MANY bins
#can change bin width 

ggplot(data = dvst) + 
  geom_bar(mapping = aes(x = GC.binned, colour = GC.binned))
#outlined

ggplot(data = dvst) + 
  geom_bar(mapping = aes(x = GC.binned, fill = GC.binned))
#filled

ggplot(data = dvst) + 
  geom_bar(mapping = aes(x = GC.binned, fill = cent))
#centromere added, fill

ggplot(data = dvst) + 
  geom_bar(mapping = aes(x = GC.binned, fill = cent), position = "fill")
#fill to compare PROPORTIONS only

ggplot(data = dvst) + 
  geom_bar(mapping = aes(x = GC.binned, fill = cent), position = "dodge")
#dodge to compare regions, set them side by side

# NEW DATASETS
#Read datasets
mtfs <- read_tsv("https://raw.githubusercontent.com/vsbuffalo/bds-files/master/chapter-08-r/motif_recombrates.txt")
rpts <- read_tsv("https://raw.githubusercontent.com/vsbuffalo/bds-files/master/chapter-08-r/motif_repeats.txt")
head(mtfs)
head(rpts)


rpts2 <- rpts %>% 
  unite(pos, chr, motif_start, sep = "-") %>% ## new function!
  select(name, pos) %>% 
  inner_join(mtfs, by = "pos")
#unite unites two columns in one dataset
#pos is a new column combining chromosome and motif_start, separated by a dash
#join joins two datasets based on a column