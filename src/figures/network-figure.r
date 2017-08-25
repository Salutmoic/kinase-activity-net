library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)
library(scales)
library(cowplot)
library(viridis)
library(classInt)
library(grid)
library(RColorBrewer)

citations <- read.delim("out/kinase-citations.tsv", as.is=TRUE)
citations <- subset(citations, count.type=="num.subs")

example_net1 <- ggdraw() + draw_image("img/figures/panels/network-egf-basic-akt.eps")
example_net2 <- ggdraw() + draw_image("img/figures/panels/network-egf-basic-erk.eps")
example_net3 <- ggdraw() + draw_image("img/figures/panels/network-egf-full-akt.eps")
example_net4 <- ggdraw() + draw_image("img/figures/panels/network-egf-full-erk.eps")
dummy_citations <- citations %>%
    ggplot(aes(x=kinase, y=num.cites.filt, colour=log10(num.cites.filt))) +
    geom_bar(stat="identity") +
    scale_colour_gradientn(expression(paste(log[10], "(kinase # citations)")),
                           colours=rev(brewer.pal(name="RdYlBu", n=5)),
                           values=rescale(c(0, 1, 2, 3, 4))) +
    theme(legend.justification="top", legend.title=element_text(size=7),
          legend.text=element_text(size=6),
          legend.key.height=unit(0.3, "cm"),
          legend.key.width=unit(0.35, "cm"),
          legend.direction="horizontal")
citations_legend <- get_legend(dummy_citations)
net_legend <- ggdraw() +
    draw_image("img/figures/panels/network-legend.eps", scale=2, clip="off")

save_plot("tmp/net-legend.pdf", net_legend)

top_grid <- plot_grid(example_net1, example_net2,
                      nrow=1, scale=0.85, labels=c("", "b"),
                      rel_widths=c(1.0, 1.0))

legends_grid <- plot_grid(citations_legend, net_legend, rel_heights=c(0.1, 1.0),
                          nrow=2)

blah_grid <- plot_grid(example_net3, scale=0.95)

bot_grid <- plot_grid(blah_grid, legends_grid,
                      nrow=1, scale=0.9, rel_widths=c(1.0, 0.25))

left_grid <- plot_grid(top_grid, bot_grid,
                        nrow=2, scale=1.0, rel_heights=c(1.0, 1.0),
                       labels=c("", "d"))

right_grid <- plot_grid(example_net4, scale=0.95)

net_figure <- plot_grid(left_grid, right_grid, nrow=1,
                        scale=1.0, rel_widths=c(1.0, 1.0),
                        labels=c("a", "c"))

save_plot("img/figures/network-figure.pdf", net_figure, ncol=2,
          nrow=1)
