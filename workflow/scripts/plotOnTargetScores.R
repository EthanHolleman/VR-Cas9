library(ggplot2)
library(ggpubr)
library(RColorBrewer)


read_table <- function(table.path){

    read.table(file = table.path, sep = '\t', header = TRUE)

}

make_plot <- function(table){


    scatter <- ggplot(table, aes(x=start, y=Doench2014OnTarget, fill=selected)) +
                geom_point(color='black', pch=21, size=3) + theme_pubr() + 
                scale_fill_brewer(palette='Dark2')
    
    box <- ggplot(table, aes(x=selected, y=Doench2014OnTarget, fill=selected)) + geom_boxplot() +
        theme_pubr() + stat_compare_means(method = "t.test", label = "p.signif") +
        scale_fill_brewer(palette='Dark2') + theme(legend.position = "none")
    
    ggarrange(scatter, box, widths = c(2, 0.8))

}

main <- function(){

    table <- read_table(snakemake@input[[1]])
    plot <- make_plot(table)
    ggsave(snakemake@output[[1]], plot, dpi=300, width=10, height=5)

}

if (! interactive()){

    main()

}

