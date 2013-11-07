$(function() {
  
  var gridster = $(".gridster ul").gridster().data('gridster');
  var gridster_enabled = true;
  var gridster_controls_show = true;
  
  $("#gridster-control").click( function() {
    
    if (gridster_enabled) {
      gridster.disable();
      gridster_enabled = false;
      $("#gridster-control").html("Gridster Disabled");
    } else {
      gridster.enable();
      gridster_enabled = true;
      $("#gridster-control").html("Gridster Enabled");
    }
    
  });
  
  $("#gridster-control-hide").click( function() {
    if (gridster_controls_show) {
      $("#controls-container").css("display", "inherit");
      $("#gridster-control-container").animate({
        height: "95%"
      });
      // $("body").css("overflow", "hidden"); // used if we want to disable scrolling
      $("#gridster-control-hide").html("Hide Controls");
      gridster_controls_show = false;
    } else {
      $("#controls-container").css("display", "none");
      $("#gridster-control-container").animate({
        height: "25px"
      });
      // $("body").css("overflow", "auto"); // used if we want to re-enable scrolling
      $("#gridster-control-hide").html("Show Controls");
      gridster_controls_show = true;
    }
  });
  
});
