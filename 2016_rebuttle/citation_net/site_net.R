library(igraph)

get_useful_col <- function(w){
  w[,c("PT","AU","TI","SO","DT","PY","DI","WC","UT")]
}


get_citation <- function(string_vec){
  first_author <- string_vec[1] |>
    stringr::str_split(",")
  first_author <- first_author[[1]][1]
  if(length(string_vec)==1) return(first_author)
  if(length(string_vec)==2){
    second_author <- string_vec[2] |>
      
      stringr::str_split(",")
    second_author <- second_author[[1]][1]
    return(paste0(first_author," &",second_author))
  }
  
  return(paste(first_author,"et al."))
}

is_self_citation <- function(from, to = "Chapron, G; Treves, A"){
  author_list_to <- stringr::str_split(to,"; ")[[1]]
  author_list_from <- stringr::str_split(from,"; ") |>
    sapply(FUN = function(w, author_list_to){
      length(intersect(w, author_list_to ))>0
    }, author_list_to)
  author_list_from
}

make_igraph <- function(all_edges, authorship){
  citation_graph <- graph_from_data_frame(d = all_edges, vertices = authorship)
  color <- c("black","red")
  E(citation_graph)$color <- color[E(citation_graph)$self_citation + 1]
  V(citation_graph)$name <- authorship[V(citation_graph)$name,"citation"]
  return(citation_graph)
}


treves <- read.table("./savedrecs.txt", header = T,sep = "\t", quote = "") |> get_useful_col()
stien <- read.table("./stien_comments.txt", header = T, sep = "\t", quote = "") |> get_useful_col()
olson <- read.table("./olson_comments.txt", header = T, sep = "\t", quote = "") |> get_useful_col()
pepin <-  read.table("./pepin_comments.txt", header = T, sep = "\t", quote = "") |> get_useful_col()

authorship <- rbind(treves[,c("UT","AU","PY","TI")],stien[,c("UT","AU","PY","TI")],olson[,c("UT","AU","PY","TI")],pepin[,c("UT","AU","PY","TI")]) |>
  unique()

citations <- authorship$AU |> stringr::str_split(";") |>
  sapply(FUN = get_citation) |>
  paste(authorship$PY)

authorship$citation <- citations
row.names(authorship) <- authorship$UT
order_papers <- c("WOS:000391104300019","WOS:000416391400008","WOS:000416391400007","WOS:000397884000001",authorship$UT) |>
  unique()
authorship <- authorship[order_papers,]


edge_treves <- data.frame(from = treves$UT[treves$UT != "WOS:000391104300019"],to = "WOS:000391104300019", self_citation = is_self_citation(from = treves$AU[treves$UT != "WOS:000391104300019"]))
edge_stien <- data.frame(from = stien$UT,to = "WOS:000416391400008", self_citation = is_self_citation(from = stien$AU, to = "Stien, A"))
edge_olson <- data.frame(from = olson$UT,to = "WOS:000416391400007", self_citation = is_self_citation(from = olson$AU, to = "Olson, ER; Crimmins, SM; Beyer, DE; MacNulty, DR; Patterson, BR; Rudolph, BA; Wydeven, AP; Van Deelen, TR"))
edge_pepin <- data.frame(from = pepin$UT,to = "WOS:000397884000001", self_citation = is_self_citation(from = pepin$AU, to = "Pepin, KM; Kay, SL; Davis, AJ"))

all_edges <- rbind(edge_treves,edge_stien,edge_olson,edge_pepin)

#g_treves <- make_igraph(edge_treves, authorship)
#g_stien <- make_igraph(edge_stien, authorship)
#g_olson <- make_igraph(edge_olson, authorship)
#g_pepin <- make_igraph(edge_pepin, authorship)

g_combine <- make_igraph(all_edges, authorship)
V(g_combine)$color <- c("orange","chartreuse")[(degree(g_combine) >= 2)+1]
V(g_combine)$color[1:4] <- c("gray",rep("lightblue",3))
V(g_combine)$color[67] <- "cyan1"

#g_combine <- igraph::union(g_treves, g_stien)
prot_star <- make_star(64, "in")
pro_coords <- layout_as_star(prot_star)
in_the_center <- matrix(c(-.3,.3,.3,.3,-.3,-.3,.3,-.3), ncol = 2,byrow = T)
coords <- rbind(in_the_center, pro_coords[-1,])


png("cite_net.png", width = 8, height = 5, res = 500, unit = "in")
par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0), xpd = TRUE)
plot(g_combine, layout = 2*coords, edge.arrow.size = 0.15, 
     edge.width= 1,vertex.size=c(rep(40,4),rep(8,63)), 
     vertex.label = c("Chapron\nand\nTreves\n2016", "Stien \n2017", "Olson et al. \n2017",  "Pepin et al. \n2017"  , 1:63),
     xlim = c(-1,1), ylim = c(-.8,.8)
     )
legend("topleft",legend = c("cites original","cites both","cites critic","self citation"),pt.bg = c("orange","chartreuse","cyan1",NA),col = c(NA,NA,NA,"red"),lty = c(0,0,0,1),pt.cex = 2,pch = c(22,22,22,NA), cex = 0.8)
dev.off()

# alt score rank
trevers_alt <- read.table("./savedrecs_alt.txt", header = T,sep = "\t", quote = "") 
rownames(trevers_alt) <- trevers_alt$UT
trevers_alt <- trevers_alt[order_papers,]
#paper_ids <- c("Chapron & Treves 2016", "Stien 2017", "Olson et al. 2017",  "Pepin et al. 2017"  , 1:63)
paper_ids <- c("C&T2016", "S2017", "O et al. 2017",  "P et al. 2017"  , 1:63)
trevers_alt$id <- paper_ids
alt_order <- rev(order(trevers_alt$Altmetric))
trevers_alt$id <- factor(trevers_alt$id, levels = paper_ids[alt_order])
trevers_alt <- trevers_alt[!is.na(trevers_alt$Altmetric), ]
trevers_alt$AU <- sub('"',"",sub('"',"",trevers_alt$AU))
colors <- c("gray","black","red")
names(colors) <- c("original","No","Yes")
self_cite <- is_self_citation(trevers_alt$AU) + 2
trevers_alt$self_cite <- c("original","No","Yes")[c(1,self_cite[-1])]
trevers_alt$self_cite <- factor(trevers_alt$self_cite, levels = c("original","No","Yes"))
library(ggplot2)
ggplot(trevers_alt, aes(x = id, y = Altmetric)) + 
  geom_bar(stat="identity",aes(fill = self_cite)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = -65, vjust = 0.5, hjust=0)) + 
  xlab("") + 
  labs(fill = "Self citation") +
  scale_fill_manual(values = colors) + 
  theme(legend.position = c(0.9, 0.7))

  
# 
ggsave("altmetric.png", width = 10, height = 6, unit = "in", dpi = 500, scale = 0.8)


