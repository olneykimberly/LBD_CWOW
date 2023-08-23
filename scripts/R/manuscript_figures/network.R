library(visNetwork)
install.packages("ggnetwork")
library(ggnetwork)
data <- toVisNetworkData(magenta_graph)
nodes <- data$nodes
edges <- data$edges

igraph::graph.adjacency(data)
ggplot(ggnetwork(magenta_graph))

visNetwork(nodes, edges, width = "100%") %>% 
  visNodes(shape = "square", 
           color = list(background = "lightblue", 
                        border = "darkblue",
                        highlight = "yellow"),
           shadow = list(enabled = TRUE, size = 10)) %>%
  visLayout(randomSeed = 12) # to have always the same network 


visNetwork(nodes, edges, height = "500px", width = "100%") # dancing plot

nodes <- data.frame(id = 1:10,
                    
                    # add labels on nodes
                    label = paste("Node", 1:10),
                    
                    # add groups on nodes 
                    group = c("GrA", "GrB"),
                    
                    # size adding value
                    value = 1:10,          
                    
                    # control shape of nodes
                    shape = c("square", "triangle", "box", "circle", "dot", "star",
                              "ellipse", "database", "text", "diamond"),
                    
                    # tooltip (html or character), when the mouse is above
                    title = paste0("<p><b>", 1:10,"</b><br>Node !</p>"),
                    
                    # color
                    color = c("darkred", "grey", "orange", "darkblue", "purple"),
                    
                    # shadow
                    shadow = c(FALSE, TRUE, FALSE, TRUE, TRUE))             

# head(nodes)
# id  label group value    shape                     title    color shadow
#  1 Node 1   GrA     1   square <p><b>1</b><br>Node !</p>  darkred  FALSE
#  2 Node 2   GrB     2 triangle <p><b>2</b><br>Node !</p>     grey   TRUE

edges <- data.frame(from = c(1,2,5,7,8,10), to = c(9,3,1,6,4,7))
visNetwork(nodes, edges, height = "500px", width = "100%")



#--------------------------------------------------------------------------------

library("grid")
library("ggplotify")
ggplot(magenta_graph, edge.arrow.size=.2, edge.color="magenta",
       vertex.color="magenta", vertex.frame.color="#ffffff",
       vertex.label=V(magenta_graph)$gene_name, vertex.label.color="black")
#path <- paste0("../../../results/manuscript_figures/magenta")
#saveToPDF(paste0(path, ".pdf"), width = 4, height = 4)
graph_p
#library(visNetwork)
#library(htmlwidgets)
#saveWidget(visIgraph(magenta_graph), file = "test.html") # dancing plot 

data <- toVisNetworkData(magenta_graph)



#visNetwork(nodes = data$nodes, edges = data$edges, height = "500px")
kitty <- visNetwork(data$nodes, data$edges, height = "500px") %>%
  visIgraphLayout() %>%
  visNodes(size = 10)

visNetwork(data$nodes, data$edges, height = "500px") %>%
  visIgraphLayout(layout = "layout_in_circle") %>%
  visNodes(size = 10) %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), 
             nodesIdSelection = T)


mygraph <- igraph::graph_from_data_frame(magenta_graph) #create a graph
mylayout <- igraph::layout_as_tree(magenta_graph) #create a circular layout

g <- make_ring(10) + make_full_graph(5)
coords <- layout_(magenta_graph)
plot(magenta_graph, layout = coords)

corrds <- layout_nicely(magenta_graph, dim = 2)

magenta_graph$layout = mylayout #store layout as a graph attribute


#another gotcha is that ig2ggplot needs both vertex names and vertex labels. 
#as of now you just have vertex names. 
V(magenta_graph)$label = V(magenta_graph)$name #store label as a vertex attrbute

#install.packages("remotes")
#remotes::install_github("etheleon/MetamapsDB")
library(MetamapsDB)

ig2ggplot(magenta_graph, dfOnly = TRUE, labels = gene_name, metab = TRUE)

ggplot(data)
typeof(magenta_graph)

MetamapsDB::ig2ggplot(DEG_graph, 
                      dfOnly = FALSE, 
                      labels = FALSE, 
                      metab = TRUE ) + 
  theme(legend.position = 'none')
mytest_plot <- MetamapsDB::ig2ggplot(magenta_graph, 
                                     dfOnly = FALSE, 
                                     labels = FALSE, 
                                     metab = TRUE ) + theme(legend.position = 'none')



test <- ggraph(DEG_graph, layout = 'kk') + 
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) + 
  geom_node_point(aes(colour = factor(selection))) + 
  theme_graph(foreground = 'white', fg_text_colour = 'white') 


lay <- layout_nicely(DEG_graph)
lay <- as.matrix(lay)

plot(DEG_graph, edge.arrow.size=.2, edge.color="magenta" ,
     vertex.color="magenta" , vertex.frame.color="#ffffff",
     vertex.label=V(DEG_graph)$gene_name, vertex.label.color="black", 
     edge.label=NULL, 
     vertex.label.font = ifelse(V(DEG_graph)$selection == "TRUE", 2, 1)) 

plot(DEG_graph, vertex.label = ifelse(V(DEG_graph)$selection == "TRUE", V(DEG_graph)$gene_name, "."), vertex.label.color="black")

ig2ggplot(DEG_graph, labels = gene_name, layout = lay)

ggnetwork(DEG_graph, layout = "target", niter = 100)
test <- plot(DEG_graph, layout = lay)
ggplot(DEG_graph, layout = lay)

ggraph(DEG_graph, layout = 'nicely',) + 
  theme_graph()


