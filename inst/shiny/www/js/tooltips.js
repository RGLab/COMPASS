$( function() {

  var helpMessages = {
    
    "markers-help": {
      title: "Select Subsets Expressing the Selected Markers",
      message: "This selection filters the subsets presented in the heatmap, " +
        "as well as available in the lower <strong>Visualize subsets...</strong> " +
        "dropdown, to those expressing <strong>all</strong> of the selected " +
        "subsets."
    },
    
    "dof-help": {
      title: "Filter Subsets by Degree of Functionality",
      message: "The degree of functionality is the number of markers expressed " +
        "in a particular subset. Susbets with a degree of functionality outside " +
        "of the specified range will not be plotted in the heatmap."
    },
    
    "facets-help": {
      title: "Conditioning Variables",
      message: "The metadata variables selected here can be used to facet plots " +
        "in different ways. This control facilitates comparison of cellular functionality " +
        "by the selected metadata variables."
    },
    
    "filter1-help": {
      title: "Filter Variables",
      message: "The data presented in the various plots can be filtered according " +
        "to different metadata information available."
    },
    
    "subsets-help": {
      title: "Select Subsets for Visualization in the Histogram",
      message: "The subsets selected here will be visualized in the histogram."
    }
    
  };
  
  $(".fa-question-circle").each( function(index, value) {
    console.log("foo!");
    var id = $(this).attr("id");
    var msg = helpMessages[id];
    var help_message = "<strong>" + msg.title + "</strong><br/><br/>" + msg.message;
    new Opentip("#" + id, help_message, { fixed: true, background: "#FFFFBA", borderColor: "#EAEAEA" });
  })

});
