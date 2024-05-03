#shiny server


server <- shinyServer(function(input, output, session) {
  
  #generate logo
  output$logo <- renderImage({
    list(src = '../../data/dna_tree_of_life.jpg', height = 50, width=70)
    
  }, deleteFile = F)
  
  #import preprocessed data
  #master_tbl <<- read_csv('/home/dylan/Postdoc/Hindle/broad/data/cleaned_EVE_data.csv')
  thetas <<- read_csv('../../data/shinyEveData/cleaned_EVE_thetas.csv')
  exprs <<- read_csv('../../data/shinyEveData/cleaned_EVE_expr.csv')
  lfc <<- read_csv('../../data/shinyEveData/cleaned_EVE_lfc.csv')
  
  gname_trans <<- thetas %>% select(gene, clean_gname) %>% unique()
  
  tree <- read.tree('../../data/shinyEveData/9_species-list.nwk')
  
  
  #populate gene_select list
  output$gene_select <- renderUI({
    selectizeInput("igene", "Select a Gene", choices = gname_trans %>% pull(clean_gname) %>% unique(), selected = 'GLT8D2')
  })
  
  #banned samples
  banned_samps <<- c('human_HS2096_37C_2-5mM',
                   'human_HS2096_37C_30mM',
                   'human_HS3099_37C_2-5mM',
                   'human_HS3099_37C_30mM',
                   'human_HS3896_37C_2-5mM',
                   'human_HS3896_37C_30mM',
                   'human_HS4057_37C_2-5mM',
                   'human_HS4057_37C_30mM',
                   'human_HS967_37C_2-5mM',
                   'human_HS967_37C_30mM')
  
  #prepare expression table 
  exprs <<- exprs %>% 
    select(-any_of(banned_samps) ) %>% 
    select(gene, contains('_8mM')) %>% 
    pivot_longer(cols = contains('_8mM'), names_to = 'sample_id', values_to = 'tpm_expr') %>% 
    mutate(species = case_when(
      grepl('human', sample_id) ~ 'human',
      grepl('mouse', sample_id) ~ 'mouse',
      grepl('bactrian_camel', sample_id) ~ 'bactrian_camel',
      grepl('dromedary_camel', sample_id) ~ 'dromedary_camel',
      grepl('bat', sample_id) ~ 'bat',
      grepl('dolphin', sample_id) ~ 'dolphin',
      grepl('honey_badger', sample_id) ~ 'honey_badger',
      grepl('rhino', sample_id) ~ 'rhino',
      grepl('rat', sample_id) ~ 'rat',
      grepl('squirrel', sample_id) ~ 'squirrel',
    )) %>%
    mutate(temp = case_when(
      grepl('32C', sample_id) ~ '32C',
      grepl('37C', sample_id) ~ '37C',
      grepl('41C', sample_id) ~ '41C'
    )) %>%
    merge(., gname_trans, by='gene')
  
  #prepare lfc table
  lfc <<- lfc %>% 
    select(-any_of(banned_samps) ) %>% 
    select(gene, contains('_8mM')) %>% 
    pivot_longer(cols = contains('_8mM'), names_to = 'sample_id', values_to = 'tpm_lfc') %>% 
    mutate(species = case_when(
      grepl('human', sample_id) ~ 'human',
      grepl('mouse', sample_id) ~ 'mouse',
      grepl('bactrian-camel', sample_id) ~ 'bactrian_camel',
      grepl('dromedary-camel', sample_id) ~ 'dromedary_camel',
      grepl('bat', sample_id) ~ 'bat',
      grepl('dolphin', sample_id) ~ 'dolphin',
      grepl('honey-badger', sample_id) ~ 'honey_badger',
      grepl('rhino', sample_id) ~ 'rhino',
      grepl('rat', sample_id) ~ 'rat',
      grepl('squirrel', sample_id) ~ 'squirrel',
    )) %>%
    mutate(temp = case_when(
      grepl('32C', sample_id) ~ '32C',
      grepl('37C', sample_id) ~ '37C',
      grepl('41C', sample_id) ~ '41C'
    )) %>%
    merge(., gname_trans, by='gene')
  
  
  #refresh res_table for output
  observeEvent(input$refresh_table, {
    req(input$igene)
    req(input$iclade)
    req(input$itemp)
    print('refresh')
    print(c(input$igene, input$iclade, input$itemp, input$idatalevel))
    print(thetas)
    output$res_table <- DT::renderDataTable({
      DT::datatable(
        thetas %>% 
          filter(clade == input$iclade) %>% 
          filter(temp == input$itemp) %>% 
          filter(datalevel==input$idatalevel) %>%
          arrange(-LRT) %>% 
          relocate(padjust, .after = LRT) %>%
          relocate(pval, .after = LRT) %>% 
          relocate(clean_gname), 
        list(pageLength = 15, scrollX=T)
      )
      
    })
    
    
  })
  
  observeEvent(input$refresh_plot, {
    req(input$igene)
    req(input$iclade)
    req(input$itemp)
    req(input$idatalevel)
    
    #render the main plot
    output$plot <- renderPlot({
      
      #first prepare expr / lfc table for plotting
      if(input$idatalevel=='exp.level'){
      tbl <- exprs %>% 
        filter(temp == input$itemp) %>% 
        filter(clean_gname == input$igene) %>%
        mutate(species = factor(species, levels = c('human', 
                                                    'squirrel',
                                                    'rat', 
                                                    'bat', 
                                                    'honey_badger',
                                                    'rhino', 
                                                    'dolphin', 
                                                    'dromedary_camel',
                                                    'bactrian_camel')),
               xstat = log2(tpm_expr),
               control=F) %>%
        unique()
      } else { #fold change
        tbl <- lfc %>% 
          filter(temp == input$itemp) %>% 
          filter(clean_gname == input$igene) %>%
          mutate(species = factor(species, levels = c('human', 
                                                      'squirrel',
                                                      'rat', 
                                                      'bat', 
                                                      'honey_badger',
                                                      'rhino', 
                                                      'dolphin', 
                                                      'dromedary_camel',
                                                      'bactrian_camel')),
                 xstat = tpm_lfc,
                 control=F) %>%
          unique()
      }
      
      #calc theta/shift values
      theta <- thetas %>% 
        filter(datalevel == input$idatalevel) %>%
        filter(clade == input$iclade) %>% 
        filter(temp == input$itemp) %>% 
        filter(clean_gname == input$igene) %>% 
        pull(theta)
      thetashift <- thetas %>% 
        filter(datalevel == input$idatalevel) %>%
        filter(clade == input$iclade) %>% 
        filter(temp == input$itemp) %>% 
        filter(clean_gname == input$igene) %>% 
        pull(thetaShift)
      if (theta<thetashift){
        vjust_theta <- -0.7
        vjust_thetashift <- 1.2
      } else {
        vjust_theta <- 1.2
        vjust_thetashift <- -0.7
      }
      
      #check plot settings 
      if ((input$show_control=='Yes') & (input$idatalevel=='exp.level')) {
        
        ##append control data 
        tmp <- exprs %>% 
          filter(temp == "37C") %>% 
          filter(clean_gname == input$igene) %>%
          mutate(species = factor(species, levels = c('human', 
                                                    'squirrel',
                                                    'rat', 
                                                    'bat', 
                                                    'honey_badger',
                                                    'rhino', 
                                                    'dolphin', 
                                                    'dromedary_camel',
                                                    'bactrian_camel')),
                xstat = log2(tpm_expr),
                control=T) %>%
          unique()

          tbl <- rbind(tbl %>% mutate(control=F), tmp)
          
          if (input$iclade=='flexible') {
            ##create custom scale fill manual
            fill_vals <- c('#42eff5', 'grey', #human
                          '#fc8c86', 'grey', #squirrel
                          '#42eff5', 'grey', #rat
                          '#fc8c86', 'grey', #bat
                          '#42eff5', 'grey', #honey badger
                          '#42eff5', 'grey', #rhino
                          '#42eff5', 'grey', #dolphin
                          '#fc8c86', 'grey', #drom camel
                          '#fc8c86', 'grey') #bac camel
          } else {
            fill_vals <- c('#42eff5', 'grey', #human
                          '#fc8c86', 'grey', #squirrel
                          '#42eff5', 'grey', #rat
                          '#fc8c86', 'grey', #bat
                          '#42eff5', 'grey', #honey badger
                          '#42eff5', 'grey', #rhino
                          '#42eff5', 'grey', #dolphin
                          '#42eff5', 'grey', #drom camel
                          '#42eff5', 'grey') #bac camel
          }
          
          point_col_vals <- c('#c96ffc', 'black')
          shape_scale_vals <- c(17, 16)
          x_lab <- paste0('Log2 TPM ' , input$itemp)
      } else { #meaning either lfc data or no control data
        if (input$iclade=='flexible'){
          fill_vals <- c('#42eff5', #human
                        '#fc8c86', #squirrel
                        '#42eff5', #rat
                        '#fc8c86', #bat
                        '#42eff5', #honey badger
                        '#42eff5', #rhino
                        '#42eff5', #dolphin
                        '#fc8c86', #drom camel
                        '#fc8c86') #bac camel
        } else {
          fill_vals <- c('#42eff5', #human
                         '#fc8c86', #squirrel
                         '#42eff5', #rat
                         '#fc8c86', #bat
                         '#42eff5', #honey badger
                         '#42eff5', #rhino
                         '#42eff5', #dolphin
                         '#42eff5', #drom camel
                         '#42eff5') #bac camel
        }
        point_col_vals <- c('black')
        shape_scale_vals <- c(16)
        
        if (input$idatalevel == 'exp.level'){
          x_lab <- paste0('Log2 TPM ' , input$itemp)
        } else {
          x_lab <- paste0('Log2 Fold Change (Gene Score) ' , input$itemp, 'vs 37C')
        }
        
      }
      
      #plot the tree
      left <- ggtree(tree) + coord_cartesian(clip="off") #+ geom_tiplab() 
      
      right <- ggplot(tbl, mapping=aes(x=xstat, y=species, fill=interaction(control,species), color=control, shape=control))+
          geom_boxplot(color = 'black',outlier.shape = NA)+
          geom_jitter(alpha=0.6, height=0.1, width=0.1)+
          scale_fill_manual(values = fill_vals)+ 
          scale_color_manual(values = point_col_vals)+
          scale_shape_manual(values = shape_scale_vals)+
          theme_minimal()+
          theme(
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_text(face='bold', size=12),
            axis.line.y=element_blank(),
          )+
          xlab(x_lab)+
          guides(fill='none', color='none', shape='none')
      

      #show theta as vlines?
      if (input$show_thetas=='Yes') {
        right <- right+
          geom_textvline(label = 'Background Theta', xintercept = theta, linetype='dashed', vjust=vjust_theta)+
          geom_textvline(label = 'Shifted Theta',xintercept = thetashift, linetype='dotted', vjust=vjust_thetashift)
        
      }
      plot_grid(left, right, ncol=2, align='h', rel_widths=c(3,9))
    })
    
    ##render a bunch of statistics plots for this gene
    output$sidebar_plot <- renderPlot({
      #filter data
      stbl <- thetas %>%
        filter(clade == input$iclade) %>% 
        filter(temp == input$itemp) %>%
        filter(datalevel==input$idatalevel) %>%
        select(clean_gname, LRT, theta, thetaShift, padjust, pval, beta) 

      #find the gene and mark it
      stbl <- stbl %>% 
        mutate(label = case_when(
          clean_gname == input$igene ~ 'yes',
          TRUE ~ 'no'
        ))
      
      lrt_bounds <- log10((stbl %>% filter(clean_gname == input$igene) %>% pull(LRT)) * c(0.90,1.05)) + c(-0.1,0.1)
      padj_bounds <- (1-log10((stbl %>% filter(clean_gname == input$igene) %>% pull(padjust))) * c(0.99,1.01)) + c(-0.05,0.05)
      beta_bounds <- log10(stbl %>% filter(clean_gname == input$igene) %>% pull(beta)) * c(0.95,1.05) + c(-0.2,0.2)
      theta_bounds <- (stbl %>% filter(clean_gname == input$igene) %>%  mutate(dtheta=theta-thetaShift)%>% pull(dtheta)) * c(0.95,1.05) + c(-0.1,0.1)
      
      g1 <- ggplot(stbl, aes(x=log10(LRT)))+
        stat_density(geom='line')+
        stat_density(geom='area', aes(x=stage(log10(LRT), after_stat = oob_censor(x,lrt_bounds))), fill='red')+
        theme_minimal()+
        theme(
          axis.text.y=element_blank()
        )+
        ylab('Density')
      
      g2 <- ggplot(stbl,aes(x=1-log10(padjust))) +
        stat_density(geom='line')+
        stat_density(geom='area', aes(x=stage(1-log10(padjust), after_stat = oob_censor(x, padj_bounds))), fill='red')+
        theme_minimal()+
        theme(
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank()
        )
      
      g3 <- ggplot(stbl,aes(x=log10(beta)))+
        stat_density(geom='line')+
        stat_density(geom='area', aes(x=stage(log10(beta), after_stat = oob_censor(x, beta_bounds))), fill='red')+
        theme_minimal()+
        theme(
          axis.text.y=element_blank()
        )+
        ylab('Density')
      
      g4 <- ggplot(stbl,aes(x=theta-thetaShift))+
        stat_density(geom='line')+
        stat_density(geom='area', aes(x=stage(theta-thetaShift, after_stat = oob_censor(x, theta_bounds))), fill='red')+
        scale_x_continuous(limits = c(-10,10), expand=c(0,0))+
        theme_minimal()+
        theme(
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank()
        )
      
      body <- plot_grid(g1,g2,g3,g4, ncol=2, nrow=2)
      title <- ggdraw() + 
        draw_label(
          paste0('Gene: ', input$igene),
          fontface = 'bold',
          x = 0,
          hjust = -2.5
        )
      plot_grid(title, body, ncol=1, rel_heights=c(0.1,0.9))
    })
  })
  
  
} )
