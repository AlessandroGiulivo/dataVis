#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel(HTML('<h1>Old Faithful Geyser Data</h1>
        <hr><br>
        <p> I`m a paragraph showing how to write
        <strong>bold</strong> and <em>italics</em>.<br>
        This is a <code>block of code</code>.<br>
        I can also put <a href="http://www.google.com">a link to google</a>.<br>
        And I can also add images!<br>
        <center><img src="https://ih0.redbubble.net/image.543360195.2115/pp,550x550.jpg" width=10%></center></p>
        <h2>The app</h2>')),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
