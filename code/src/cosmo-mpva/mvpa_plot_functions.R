#plotting functions for MVPAs 

unimodal_mvpa_plot_ReadSpeech <- function(plot_data, summary_data){
    
    #Extract mean and se accuracy from summary stats
    mean_accu <- round(summary_data$mean_accuracy, digits = 2)
    se_accu <- round(summary_data$se_accuracy, digits = 2)
    
    plot_data %>% 
        ggplot(aes(x = conditions, y = accuracy, color = roiArea)) +
        geom_point(aes(alpha = I(0.7)), 
                   position = position_jitterdodge(jitter.width=0.2,
                                                   dodge.width = 0.3),
                   size = 4.5) +
        stat_summary(aes(x =conditions, y = accuracy, color = roiArea),
                     fun = "mean", 
                     geom = "crossbar", 
                     position = position_dodge(width = 0.5),
                     width = .55,
                     colour  = 'black',
                     inherit.aes = FALSE,
                     show.legend = FALSE) + 
        stat_summary(aes(x = conditions, y = accuracy, color = roiArea),
                     fun.max = function(x) mean(x) + (std.error(x)), 
                     fun.min = function(x) mean(x) - (std.error(x)), 
                     geom = "errorbar",
                     position = position_dodge(width = 0.5),
                     width = .15,
                     size = 0.8,
                     colour  = 'black',
                     inherit.aes = FALSE,
                     show.legend = FALSE) + 
        geom_hline(yintercept=c(50), 
                   linetype="dotted", 
                   colour="black", 
                   linewidth=.5) +
        #scale_color_manual(values = c(rep(dark_colors[1], 3))) +
        scale_y_continuous(name = "Clasisfier Accuracy [%]", 
                           limits = c(0,100)) +
        scale_x_discrete(name = 'Decoding conditions') +
        annotate("text", 
                 x = classifier, 
                 y = 5, 
                 label = mean_accu)+
        annotate("text", 
                 x = classifier, 
                 y = 0, 
                 label = se_accu) +
        theme_cowplot(font_size = 16, font_family = "Arial") +
        theme(axis.line = element_line(colour = 'black', size = 1),
              axis.ticks = element_line(colour = 'black', size = 1),
              axis.text = element_text(face="bold"),
              legend.position = 'none')
    
}

crossmodal_mvpa_plot_ReadSpeech <- function(plot_data, summary_data){
    
    #Extract mean and se accuracy from summary stats
    mean_accu <- round(summary_data$mean_accuracy, digits = 2)
    se_accu <- round(summary_data$se_accuracy, digits = 2)
    
    plot_data %>% 
        ggplot(aes(x = conditions, 
                   y = accuracy, 
                   color = TrainTest, 
                   alpha = TrainTest,
                   size = TrainTest,
                   group = interaction(conditions, TrainTest))) +
        geom_point(position = position_jitterdodge(jitter.width=0.2,
                                                   dodge.width = 0.9)) +
        
        #scale_alpha_manual(values = c(0.2,0.8, 0.2), 
        scale_alpha_manual(values = c(0.95,0.95, 0.95), 
                           name = "Decoding partitions") +
        # scale_color_manual(values = c(rep(dark_colors[1], 3)),
        #scale_color_manual(values = blind_triple_colors,
        #                   name = "Decoding partitions")+
        scale_size_manual(values = c(1.75,4,1.75),
                          name = "Decoding partitions")+
        scale_y_continuous(name = "Clasisfier Accuracy [%]", 
                           limits = c(0,100)) +
        scale_x_discrete(name = 'Decoding conditions') +
        ## ADD STAT SUMMARIES 
        stat_summary(fun = "mean",
                     geom = "crossbar",
                     position = position_dodge(width = 0.9),
                     width = .75,
                     size = 0.75,
                     #color = 'black',
                     alpha = 1,
                     show.legend = FALSE) +
        stat_summary(fun.max = function(x) mean(x) + (sd(x)/sqrt(20)),
                     fun.min = function(x) mean(x) - (sd(x)/sqrt(20)),
                     geom = "errorbar",
                     position = position_dodge(width = 0.9),
                     width = .15,
                     size = 0.8,
                     alpha = 1,
                     color = "black",
                     show.legend = FALSE) +
        #ADD CHANCE LINE
        geom_hline(yintercept=c(50), 
                   linetype="dotted", 
                   colour="black", 
                   linewidth=.5) +
        # ## TRY TO ANNOTATE THE PLOTS WITH MEAN AND SEM VALUES - only average!? 
        annotate("text",
                 x = 1:3,
                 y = 5,
                 label = mean_accu) +
        annotate("text",
                 x = 1:3,
                 y = 0,
                 label = se_accu) +
        theme_cowplot(font_size = 16, font_family = "Arial") +
        theme(axis.line = element_line(colour = 'black', size = 1),
              axis.ticks = element_line(colour = 'black', size = 1),
              axis.text = element_text(face="bold"),
              legend.position="bottom",
              legend.direction = "vertical",
              legend.justification = 0.60)
    
}