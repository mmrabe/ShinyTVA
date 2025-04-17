

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)
library(pracma)
library(scales)
library(rhandsontable)
library(cowplot)
library(Rmpfr)


logsumexp <- function(logs, na.rm = FALSE) {
  if(na.rm) {
    logs <- logs[!is.na(logs)]
  } else if(any(is.na(logs))) {
    return(NA_real_)
  }
  logs <- logs[logs!=-Inf]
  if(length(logs)==0) return(-Inf)
  if(length(logs)==1) return(logs)
  if(any(is.nan(logs))) return(NaN)
  if(any(logs == Inf)) return(Inf)
  maxlog <- which.max(logs)
  logs[maxlog] + log1p(sum(exp(logs[-maxlog]-logs[maxlog])))
}

logdiffexp <- function(log_a, log_b) log_b+log1mexp(log_b-log_a)


ptvaw <- function(R, tau, K, v) {
  nR <- sum(R)
  SmR <- which(!R)
  if(K < nR) {
    -Inf
  } else if(tau <= 0 || K == 0) {
    if(nR == 0) 0 else -Inf
  } else if(nR < K) {
    sum(pexp(tau, v[R], TRUE, TRUE)) + sum(pexp(tau, v[!R], FALSE, TRUE))
  } else {
    logsumexp(vapply(which(R), \(i) {
      Rmi <- R
      Rmi[i] <- FALSE
      Rmi <- which(Rmi)
      nRmi <- length(Rmi)
      sp <- double(0)
      sm <- double(0)
      for(k in 0:nRmi) {
        PRmi <- if(k > 0) combinations(nRmi, k, Rmi) else matrix(integer(), 1, 0)
        for(j in seq_len(nrow(PRmi))) {
          J <- PRmi[j,]
          vsum <- (v[i] + sum(v[J]) + sum(v[SmR]))
          if(length(J) %% 2L) {
            sp <- c(sp, pexp(tau,vsum,log.p=TRUE)-log(vsum))
          } else {
            sm <- c(sm, pexp(tau,vsum,log.p=TRUE)-log(vsum))
          }
        }
      }
      log(v[i]) + logdiffexp(logsumexp(sp), logsumexp(sm))
    }, double(1)))
  }
}


ptvap <- function(RT, tau, SD, K, v) {
  if(tau <= 0) {
    return(if(any(RT)) -Inf else 0)
  }
  nRT <- sum(RT)
  if(is.logical(SD)) {
    SD <- which(SD)
  }
  nSD <- length(SD)
  if(nSD == 0) {
    ptvaw(RT, tau, K, v)
  } else if(nRT <= K) {
    logsumexp(vapply(0:pmin(K-nRT,nSD), \(k) {
      RDs <- if(k > 0) combinations(nSD, k, SD) else matrix(integer(), nrow = 1)
      logsumexp(vapply(seq_len(nrow(RDs)), \(i) {
        R <- RT
        R[RDs[i,]] <- TRUE
        ptvaw(R, tau, K, v)
      }, double(1)))
    }, double(1)))
  } else {
    -Inf
  }
}

lscore <- function(score, exposure_duration, nS, nD, C, alpha, K, t0) {
  w <- c(rep(1,nS-nD),rep(alpha,nD))
  lchoose(nS-nD,score)+ptvap(
    c(rep(TRUE,score),rep(FALSE,nS-score)),
    exposure_duration - t0,
    c(rep(FALSE,nS-nD),rep(TRUE,nD)),
    K,
    C/1000*w/sum(w)
  )
}

lscorev <- Vectorize(lscore, c("score","exposure_duration"))

pscore <- function(...) exp(lscore(...))
pscorev <- function(...) exp(lscorev(...))

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("ShinyTVA: Shiny app for Theory of Visual Attention (TVA)"),
  
  sidebarLayout(
    sidebarPanel(
      
      inputPanel(
        tags$h4("Display settings"),
        numericInput("items",
                     "Number of items:",
                     min = 2,
                     max = 8,
                     step = 2,
                     value = 6,
                     width = "100%"),
        numericInput("distractors",
                     "Number of distractors:",
                     min = 0,
                     max = 7,
                     value = 0,
                     width = "100%"),
        numericInput("exposure_duration", "Exposure duration", min = 10, max = 500, value = 100)
      ),
      inputPanel(
        tags$h4("Model parameters"),
        sliderInput("param_C",
                    "Processing speed (\U1D436)",
                    min = 1,
                    max = 200,
                    value = 50
        ),
        sliderInput("param_t0","Sensory threshold (\U1D461\U2080)",min=-50,max=50,value=10),
        sliderInput("param_K","VSTM capacity (\U1D43E)",min=0,max=6,value=4),
         sliderInput("param_alpha",
                     "Selectivity (\U1D6FC)",
                     min = 0,
                     max = 2,
                     step = 0.1,
                     value = 0.6
         )
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Display",
          tags$h2("Display"),
          tags$p("Based on the display settings on the left, a trial could look like this:"),
          tags$hr(),
          plotOutput("display", height = "4cm", width = "4cm"),
          tags$hr(),
          tags$p(htmlOutput("displayDescription", inline = TRUE), "To randomly generate a new display using the same settings, ", actionLink("generateDisplay", "click here"),"."),
          tags$p("In a TVA experiment, display types, especially the number of distractors, will alternate between trials. Note that this display will be shown for an extremely short duration. This is often enough to report about 2-4 items. On the next tabs, you will learn how TVA explains and predicts this at different scales.")
        ),
        tabPanel(
          "Processing",
          tags$h2("Theoretical processing times"),
          tags$p("The processing capacity \U1D436 can be interpreted as the theoretical number of items that could be processed per second if memory were unlimited. Depending on the specific model implementation, TVA can make different assumptions about the distribution of \U1D436 across items. In this app, we assume it is equally distributed across locations but that targets and distractors differ with regard to their relative share, such that a distractor receives \U1D6FC-times the processing of a target. This is called ",tags$strong("filtering"),". The individual processing rates are given as: $$v_i=C\\frac{w_i}{\\sum_{j} w_j}; w_i=\\begin{cases} \\alpha & \\text{if } i \\in D  \\\\ 1 & \\text{if } i \\notin D \\end{cases} $$" %>% withMathJax()),
          tags$p("TVA assumes that all items are processed ",tags$strong("in parallel")," and that each item’s processing time is exponentially distributed, so that the probability that item \U1D456 is processed sometime before \U1D461 is given as $$\\Pr\\left(i,T\\leq t\\right)=1-\\exp(-\\tau v_i),$$ where \U1D70F is the effective exposure duration, $$\\tau=\\begin{cases} 0 & \\text{if } t < t_0 \\\\ t-t_0 & \\text{if } t \\geq t_0 \\end{cases}.$$" %>% withMathJax()),
          conditionalPanel("input.distractors > 0", tags$p(HTML("The theoretical processing time distributions for a single <font color=red>target</font> vs. a single <font color=blue>distractor</font> look as follows:"))),
          conditionalPanel("input.distractors == 0", tags$p(HTML("The theoretical processing time distribution for a single <font color=red>target</font> looks as follows:"))),
          plotOutput("singleproc"),
          tags$p("The cumulative probability formally describes the probability that the processing time for item \U1D456 is less than or equal to \U1D461. It is the same as the probability that \U1D456 was processed some time before \U1D461. Its derivative, the probability density, is more informative to compare the relative likelihood of different possible processing times. It is highest for the most likely processing time. In both panels, the vertical dashed line represents the exposure duration from the display settings pane. At that time, the model predicts ",htmlOutput("procpred",inline = TRUE)),
        ),
        tabPanel(
          "Simulate trial",
          tags$h2("Trial simulation"),
          tags$p("Based on the example display and the model parameters on the left, according to TVA, this is what could potentially happen in a single trial, which would ultimately give rise to an observable report:"),
          plotOutput("simulation", height = "3cm"),
          htmlOutput("simulationDescription", inline = FALSE),
          tags$p("To simulate another run using the parameters on the left for the same example display, ", actionLink("newSimulation", "click here"),".")
        ),
        tabPanel(
          "Scores",
          tags$h2("Scores"),
          tags$p(HTML("The model is stochastic (non-deterministic). For the same display configuration and parameter set, it can generate different report sets. Let the number of correctly reported targets be the <strong>score</strong>. Then the score for a single trial can be as low as 0 if no target was reported but it may not exceed \U1D43E or the number of targets.")),
          tags$h2("Probability distribution of scores"),
          tags$p("Based on the display settings and model parameters on the left, the theoretical distribution of predicted scores looks as follows:"),
          plotOutput("theoreticalprobs"),
          tags$p("The colored lines are the probabilities for a given score as a function of exposure duration. The solid black line is the expected (mean) score, i.e. the mean score we would expect when collecting an infinite number of trials. The vertical dashed line represents the exposure duration from the display settings pane. At that time, the model predicts ",htmlOutput("scorepred2", inline=TRUE)),
          #tags$ol(
          #  tags$li("Before \U1D461\U2080, the only possible score is 0."),
          #  tags$li("With longer exposure, the probability for an empty report gradually declines, as long as at least one target is present."),
          #  tags$li("The peaks of the score probabilities are in order, i.e. the probability for score 1 is highest before score 2 is highest etc."),
          #  tags$li("The probability of the maximum score has no visible peak; it approaches but never reaches an asymptote.")
          #)
        ),
        tabPanel(
          "Simulate subject",
          tags$h2("Simulate a subject data set"),
          tags$p("Here, you can simulate many trials with different display settings. The parameters on the left are assumed to be the same for all trials. By default, the table below lists the configuration of a typical CombiTVA paradigm with 6 items and varying numbers of distractors and exposure durations. You can edit, add, and remove settings (by right click)."),
          rHandsontableOutput("simtab"),
          tags$div(
            actionButton("simulate","Simulate"),
            downloadButton("downloadSummary","Summary")
          ),
          tags$p(HTML("By clicking on <strong>Simulate</strong>, you can run the simulation as configured above. The columns <strong>MeanScore</strong> and <strong>SEScore</strong> will be filled in with the mean and standard error of the simulated trials, respectively. You can, however, also fill in or edit those columns by hand. When the table above changes, the values are also visualized in the plot below:")),
          plotOutput("simscores"),
          tags$p("Solid lines connect mean scores for the same display type over its different exposure durations, if any. Error bars represent simple standard errors around those means. The dashed lines represent the theoretical (predicted) mean scores for the different display types. You should observe that the predicted curves align with the observed scores.")
        )
      )
    )
  )
)

generate_display <- function(nS, nD) {
  deg_per_item <- 2 * pi / nS
  deg_offset <- (nS-2)*(-0.25)*deg_per_item
  tibble(
    letter = sample(LETTERS, nS),
    target = sample(rep(c(TRUE,FALSE),c(nS-nD,nD))),
    x = cos(deg_offset+deg_per_item*0:(nS-1L)),
    y = sin(deg_offset+deg_per_item*0:(nS-1L))
  ) %>% mutate(
    color = if_else(target, "red", "blue")
  )
}

simulate_trial <- function(items, ED, C, alpha, K, t0) {
  w <- if_else(items$target, 1, alpha)
  v <- C/1000*w/sum(w)
  items %>% mutate(exposure_duration = ED, processing_time = rexp(length(v), v) + t0, processing_order = rank(processing_time), processed = processing_order <= K & processing_time <= exposure_duration, reported = processed & target)
}

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  session$userData$simulations <- data.frame(Items = rep(6L,9), Distractors = c(rep(0L,6),2:4), Exposure = as.integer(c(10,20,50,100,150,200,150,150,150)), Replications = 27L, MeanScore = rep(NA_real_,9), SEScore = rep(NA_real_,9))
  
  observe({
    updateNumericInput(inputId = "distractors", max = input$items - 1L)
    if(input$items <= input$distractors) updateNumericInput(inputId = "distractors", value = input$items - 1L)
    updateSliderInput(inputId ="param_K", max = input$items)
    if(input$param_K > input$items) updateSliderInput(inputId ="param_K", value = input$items)
  })
  
  max_score <- reactive({
    min(input$param_K, input$items-input$distractors)
  })
  
  items <- reactive({
    input$generateDisplay
    generate_display(input$items, input$distractors)
  })
  
  observeEvent(input$simtab, {
    session$userData$simulations <- hot_to_r(input$simtab)
  })
  
  observeEvent(input$simulate, {
    for(i in seq_len(nrow(session$userData$simulations))) {
      display <- generate_display(session$userData$simulations$Items[i],session$userData$simulations$Distractors[i])
      scores <- vapply(seq_len(session$userData$simulations$Replications[i]), function(j) {
        trial <- display %>% simulate_trial(session$userData$simulations$Exposure[i], input$param_C, input$param_alpha, input$param_K, input$param_t0)
        sum(trial$reported)
      }, integer(1))
      session$userData$simulations$MeanScore[i] <- mean(scores)
      session$userData$simulations$SEScore[i] <- sd(scores)/sqrt(session$userData$simulations$Replications[i])
    }
  })
  
  output$downloadSummary <- downloadHandler("Summary.csv", function(file) {
    write.csv(session$userData$simulations, file, row.names = FALSE, col.names = TRUE)
  })
  
  output$simtab <- renderRHandsontable({
    input$simulate
    rhandsontable(session$userData$simulations)
  })
  
  output$simscores <- renderPlot({
    dat <- hot_to_r(input$simtab)
    if(is.null(dat)) return(NULL)
    dat <- dat %>%
      mutate(`Display Type` = sprintf("%dT%dD", Items-Distractors,Distractors))
    predicted <- dat %>% 
      select(Items, Distractors, `Display Type`) %>% 
      unique() %>% 
      group_by_all() %>%
      rowwise() %>%
      reframe(
        Exposure = seq(0,max(dat$Exposure),length.out = 200),
        MeanScore = rowSums(as.matrix(vapply(0:min(Items,input$param_K), function(i) if(i==0) rep(0,length(Exposure)) else i*pscorev(i,Exposure,Items,Distractors,input$param_C,input$param_alpha,input$param_K,input$param_t0), double(length(Exposure)))))
      )
    
    dat %>%
      na.omit() %>%
      ggplot(aes(x=Exposure, y=MeanScore, group = `Display Type`, color = `Display Type`)) +
      theme_minimal() +
      theme(legend.position = "bottom") +
      geom_point() +
      geom_linerange(aes(ymin=MeanScore-SEScore,ymax=MeanScore+SEScore)) +
      geom_line() +
      geom_line(linetype="dashed", data=predicted) +
      labs(x = "Exposure duration", y = "Mean score", color = "Display type")
  })
  
  
  output$displayDescription <- renderText({
    targets <- items()$letter[items()$target]
    distractors <- items()$letter[!items()$target]
    if(input$distractors > 0) {
      sprintf("This is a <strong>partial-report</strong> display with %1$d target(s) (%2$s) and %3$d distractor(s) (%4$s), also referred to as <strong>%1$dT%3$dD</strong>.", input$items-input$distractors, paste("<font color=red>",targets,"</font>", sep = "", collapse=", "), input$distractors, paste("<font color=blue>",distractors,"</font>", sep = "", collapse=", "))
    } else {
      sprintf("This is a <strong>whole-report</strong> display with %1$d targets (%2$s), also referred to as <strong>%1$dT0D</strong>.", input$items, paste("<font color=red>",targets,"</font>", sep = "", collapse=", "))
    }
  })
  
  upper_bound <- reactive(max(input$exposure_duration,qexp(.95,input$param_C/1000/input$items)))
  
  single_run <- reactive({
    input$newSimulation
    simulate_trial(items(), input$exposure_duration, input$param_C, input$param_alpha, input$param_K, input$param_t0)
  })
  
  output$procpred <- renderText({
    v_target <- input$param_C/1000/(input$items+input$distractors*(input$param_alpha-1))
    v_distractor <- input$param_C/1000*input$param_alpha/(input$items+input$distractors*(input$param_alpha-1))
    if(input$distractors > 0) {
      sprintf("any given <font color=red>target</font> has been processed with a probability of %.1f%% and any given <font color=blue>distractor</font> with a probability of %.1f%%.", pexp(input$exposure_duration-input$param_t0,v_target)*100, pexp(input$exposure_duration-input$param_t0,v_distractor)*100)
    } else {
      sprintf("any given stimulus (all of which are <font color=red>targets</font>) has been processed with a probability of %.1f%%.", pexp(input$exposure_duration-input$param_t0,v_target)*100)
    }
  })
  
  output$display <- renderPlot({
    ggplot() +
      theme_void() +
      theme(plot.background = element_rect(fill = "black")) +
      annotate("point", x = 0, y = 0, color = "red", shape = 3) +
      annotate("text", x=items()$x,y=items()$y,label=items()$letter,color=items()$color,fontface="bold",size=10) +
      coord_fixed(xlim = c(-1.5,1.5), ylim = c(-1.5,1.5), expand = FALSE)
  })
  
  output$simulationDescription <- renderText({
    l <- label_ordinal()
    dat <- single_run() %>% arrange(processing_order)
    if(input$param_t0 >= input$exposure_duration) {
      sprintf("As the exposure ends (at %.0f ms after display onset) before the sensory threshold (\U1D461\U2080=%.0f ms), nothing would be reported, even if memory were available.", input$exposure_duration, input$param_t0)
    } else if(input$param_K == 0L && sum(dat$processed) == 0L) {
      "As there are no memory slots available at all (\U1D43E=0), nothing would be reported, even if items were processed during exposure."
    } else if(input$param_K == 0L && sum(dat$processed) == 1L) {
      sprintf("As there are no memory slots available at all (\U1D43E=0), nothing will be reported, even though one letter is processed during exposure.", input$param_t0)
    } else if(input$param_K == 0L && sum(dat$processed) == 2L) {
      sprintf("As there are no memory slots available at all (\U1D43E=0), nothing will be reported, even though %d letters are processed during exposure.", input$param_t0, sum(dat$processed))
    } else if(sum(dat$processed) == 0L) {
      sprintf("Processing starts at %.0f ms (\U1D461\U2080) after display onset. Even though there are memory slots available (\U1D43E=%d), no stimuli finish processing during exposure, which is why nothing will be reported.", input$param_t0, input$param_K)
    } else if(sum(dat$processed) > 0L) {
      target_msg <- if(sum(dat$reported) == 0L) {
        "No target was processed during exposure."
      } else if(sum(dat$reported) == 1L) {
        sprintf("A single target (e.g., <font color=red>%s</font>) was processed at %.0f ms.", dat$letter[dat$reported], dat$processing_time[dat$reported])
      } else {
        sprintf("There were %d targets (e.g., %s) processed during exposure, namely at %s, and they were admitted to memory %s, respectively.", sum(dat$reported), paste("<font color=red>", dat$letter[dat$reported], "</font>", sep="", collapse=", "), paste(sprintf("%.0f ms",dat$processing_time[dat$reported]), sep="", collapse=", "), paste(l(dat$processing_order[dat$reported]), sep="", collapse=", "))
      }
      distractor_msg <- if(input$distractors == 0L) {
        "Since there were no distractors in the display, none could have been processed."
      } else if(sum(dat$processed&!dat$target) == 0L) {
        "Even though there were distractors in the display, none were processed during exposure."
      } else if(sum(dat$processed&!dat$target) == 1L) {
        sprintf("A single distractor (e.g., <font color=blue>%s</font>) was processed at %.0f ms.", dat$letter[dat$processed&!dat$target], dat$processing_time[dat$processed&!dat$target])
      } else {
        sprintf("There were %d distractors (e.g., %s) processed during exposure, namely at %s, and they were admitted to memory %s, respectively.", sum(dat$processed&!dat$target), paste("<font color=blue>", dat$letter[dat$processed&!dat$target], "</font>", sep="", collapse=", "), paste(sprintf("%.0f ms",dat$processing_time[dat$processed&!dat$target]), sep="", collapse=", "), paste(l(dat$processing_order[dat$processed&!dat$target]), sep="", collapse=", "))
      }
      report_msg <- if(sum(dat$reported) == 0) {
        "Therefore no letters are reported."
      } else if(sum(dat$reported) == 1) {
        sprintf("Therefore only one letter (e.g., <strong>%s</strong>) is reported.", dat$letter[dat$reported])
      } else {
        sprintf("Therefore %d letters (e.g., <strong>%s</strong>) are reported.", sum(dat$reported), paste(sample(dat$letter[dat$reported]),collapse=""))
      }
      if(sum(dat$processed) == input$param_K) {
        if(input$param_K == 1L) {
          sprintf("Processing starts at %.0f ms (\U1D461\U2080) after display onset. The single (\U1D43E) memory slot is filled up during exposure. %s %s %s", input$param_t0, target_msg, distractor_msg, report_msg)
        } else {
          sprintf("Processing starts at %.0f ms (\U1D461\U2080) after display onset. All %d (\U1D43E) memory slots are filled up during exposure. %s %s %s", input$param_t0, input$param_K, target_msg, distractor_msg, report_msg)
        }
      } else if(sum(dat$processed) == 1L) {
        sprintf("Processing starts at %.0f ms (\U1D461\U2080) after display onset. Even though there are memory slots available (\U1D43E=%d), only one stimulus is processed during exposure. %s %s %s", input$param_t0, input$param_K, target_msg, distractor_msg, report_msg)
      } else {
        sprintf("Processing starts at %.0f ms (\U1D461\U2080) after display onset. Even though there are memory slots available (\U1D43E=%d), only %d stimuli are processed during exposure. %s %s %s", input$param_t0, input$param_K, sum(dat$processed), target_msg, distractor_msg, report_msg)
      }
    } else {
      "<font color=red>UNKNOWN SCENARIO</font>"
    }
  })
  
  output$singleproc <- renderPlot({
    wsum <- (input$param_alpha-1)*input$distractors+input$items
    plot_grid(
      ggplot() +
        theme_minimal() +
        geom_function(fun = function(x) if_else(x <= input$param_t0, 0, pexp(x-input$param_t0, input$param_C/1000/wsum)), xlim = c(input$param_t0,upper_bound()), color = "red") +
        (if(input$distractors > 0L) geom_function(fun = function(x) if_else(x <= input$param_t0, 0, pexp(x-input$param_t0, input$param_C/1000*input$param_alpha/wsum)), xlim = c(input$param_t0,upper_bound()), color = "blue", linetype = if(input$param_alpha!=1) "solid" else "dashed")) +
        annotate("segment", y = 0, x = 0, xend = input$param_t0, color = "red") +
        (if(input$distractors > 0L) annotate("segment", y = 0, x = 0, xend = input$param_t0, color = "blue", linetype = "dashed")) +
        geom_vline(xintercept = input$exposure_duration, linetype = "dashed")+
        scale_x_continuous(limits = c(0,upper_bound()), name = "Exposure duration") +
        scale_y_continuous(name = "Cumulative processing probability", limits = c(0,1)),
      ggplot() +
        theme_minimal() +
        geom_function(fun = function(x) dexp(x-input$param_t0, input$param_C/1000/wsum), xlim = c(input$param_t0,upper_bound()), color = "red") +
        (if(input$distractors > 0L) geom_function(fun = function(x) dexp(x-input$param_t0, input$param_C/1000*input$param_alpha/wsum), xlim = c(input$param_t0,upper_bound()), color = "blue", linetype = if(input$param_alpha!=1) "solid" else "dashed")) +
        annotate("segment", y = 0, x = 0, xend = input$param_t0, color = "red") +
        (if(input$distractors > 0L) annotate("segment", y = 0, x = 0, xend = input$param_t0, color = "blue", linetype = "dashed")) +
        geom_vline(xintercept = input$exposure_duration, linetype = "dashed")+
        scale_x_continuous(limits = c(0,upper_bound()), name = "Exposure duration") +
        scale_y_continuous(name = "Processing probability density"),
      nrow = 1
    )
  })
  
  output$simulation <- renderPlot({
    
    dat <- single_run()
    
    memory_capacity <- input$param_K
    exposure_duration <- input$exposure_duration
    
    
    which_first_after_ed <- if(any(dat$processing_time > exposure_duration)) min(dat$processing_order[dat$processing_time > exposure_duration]) else NA_integer_
    
    
    dat$pos <- 1+dat$processing_order+(dat$processing_order>memory_capacity)+(dat$processing_time>exposure_duration)
    K_pos <- dat$pos[match(memory_capacity,dat$processing_order)]+1
    ED_pos <- if(is.na(which_first_after_ed)) nrow(dat) + 3L else dat$pos[match(which_first_after_ed,dat$processing_order)]-1
    
    ggplot() +
      theme_void() +
      scale_x_continuous(limits=c(0,nrow(dat)+3)) +
      annotate("point",x=0,y=0) +
      annotate("point",x=dat$pos, y=0, color = if_else(dat$target, "red","blue"),shape=if_else(dat$reported, 1, 4)) +
      annotate("segment",x=1, y=0, yend=0.7, color="purple", arrow=arrow(ends = "first", length = unit(0.1, "inches"))) +
      annotate("segment",x=K_pos, y=0, yend=0.7, color="purple", arrow=arrow(ends = "first", length = unit(0.1, "inches"))) +
      annotate("segment",x=ED_pos, y=0, yend=0.7, color="purple", arrow=arrow(ends = "first", length = unit(0.1, "inches"))) +
      annotate("text",x=dat$pos,y=-0.2, label=sprintf("%.0f ms",dat$processing_time), size = 3)+
      annotate("text",x=ED_pos, y=-0.2, label=sprintf("%.0f ms",exposure_duration), size = 3, color = "purple")+
      annotate("text",x=1, y=-0.2, label=sprintf("%.0f ms",input$param_t0), size = 3, color = "purple")+
      annotate("text",x=0, y=-0.2, label=sprintf("%.0f ms",0), size = 3)+
      annotate("text",x=ED_pos, y=1, label="ED", size = 3, color = "purple")+
      annotate("text",x=K_pos, y=1, label="italic(K)", parse=TRUE, size = 3, color = "purple")+
      annotate("text",x=1, y=1, label="italic(t)[0]", parse =TRUE, size = 3, color = "purple")
    
  })
  
  
  output$theoreticalprobs <- renderPlot({
    x <- seq(input$param_t0,upper_bound(),length.out=200)
    scores <- crossing(exposure_duration = x, score = 0:max_score()) %>% mutate(p = pscorev(score,exposure_duration,input$items,input$distractors,input$param_C,input$param_alpha,input$param_K,input$param_t0)) %>% bind_rows(tibble(exposure_duration=0,score=0:max_score())%>%mutate(p=if_else(score==0,1,0)))
    expected_score <- rowSums(as.matrix(vapply(0:max_score(), function(i) if(i==0) rep(0,length(x)) else i*pscorev(i,x,input$items,input$distractors,input$param_C,input$param_alpha,input$param_K,input$param_t0), double(length(x)))))
    ggplot() +
      theme_minimal() +
      theme(legend.position = "bottom") +
      labs(x="Exposure duration", color = "Score") +
      geom_vline(xintercept = input$exposure_duration, linetype = "dashed") +
      geom_line(aes(x=exposure_duration,y=p,color=as.factor(score)), scores) +
      annotate("line",x=c(0,x),y=c(0,if(max_score() > 0) expected_score/max_score() else expected_score),color="black",linewidth=1) +
      scale_y_continuous(name = "Score probability", sec.axis = sec_axis(~.*max(1,max_score()), name = "Expected (mean) score"), limits = c(0,1))
  })
  
  output$scorepred2 <- renderText({
    p <- pscorev(0:max_score(),input$exposure_duration,input$items,input$distractors,input$param_C,input$param_alpha,input$param_K,input$param_t0)
    paste0(paste(sprintf("a score of %d with %.1f%%", 0:max_score(), p*100), collapse=", "),sprintf(", and therefore an expected score of %.1f.", sum((0:max_score())*p)))
  })
  
  
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
