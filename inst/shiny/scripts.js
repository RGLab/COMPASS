$( function() {
  
  var zoomStatsClicked = false;
  
  // we need to make sure the #stats box expands if the user clicks the expand button
  $("#zoom-stats").click( function() {
    
    var $stats = $("#stats");
    
    // get the parent height
    var height = $stats.parent().height();
    var $dataTable = $("#stats").find(".dataTable").dataTable();
    var $dtSettings = $dataTable.fnSettings();
    
    if (zoomStatsClicked) {
      
      $dtSettings._iDisplayLength = 3;
      $dataTable.fnDraw();
      
      $stats.attr("style", "height: 300px !important");
      zoomStatsClicked = false;
    } else {
      
      // guess the number of rows to display based on the height
      // TODO: no magic numbers. 37 corresponds to the row height;
      // 100 is the vertical spacing around the table
      $dtSettings._iDisplayLength = Math.max(1, Math.round( (height - 200) / 40 ));
      $dataTable.fnDraw();
      
      $stats.attr("style", "height: " + (height - 40) + "px !important");
      zoomStatsClicked = true;
    }
    
  });
  
  // multiselect
  $("#markers").multiselect({
    noneSelectedText: "Subsets must express...", 
    header: "Select markers",
    selectedText: "# of # markers selected",
    show: ['slide', 200]  
  });
  
  $("#subsets").multiselect({
    noneSelectedText: "Subsets to Visualize in Histogram...",
    header: "Select Subset",
    selectedText: "# of # subsets selected",
    show: ['slide', 200]
  });
  
  /*
  // disable body scrolling when inside the multiselect
  $(".ui-multiselect-menu").mouseover( function() {
    $("body").css("overflow", "hidden");
  });
  
  $(".ui-multiselect-menu").mouseout( function() {
    $("body").css("overflow", "initial");
  });
  */
  
  // similarly for #stats
  /*
  $("#stats").mouseover( function() {
    $("body").css("overflow", "hidden");
  });
  
  $("#stats").mouseout( function() {
    $("body").css("overflow", "initial");
  });
  */
  
  
  
});
