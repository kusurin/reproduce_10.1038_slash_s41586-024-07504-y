data <- read.csv("overview_sankey.CSV", header = FALSE)

nodes <- data.frame(name = unique(c(
  paste(data$V1, "L3", sep = "_"),
  paste(data$V2, "L2", sep = "_"),
  paste(data$V3, "L1", sep = "_")
)))

links <- data.frame(
  target = match(paste(data$V1, "L3", sep = "_"), nodes$name) - 1,
   source= match(paste(data$V2, "L2", sep = "_"), nodes$name) - 1,
  value = 1
)

links <- rbind(links, data.frame(
  target = match(paste(data$V2, "L2", sep = "_"), nodes$name) - 1,
   source= match(paste(data$V3, "L1", sep = "_"), nodes$name) - 1,
  value = 1
))

networkD3::sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target", Value = "value", NodeID = "name", units = "TWh", fontSize = 12)
