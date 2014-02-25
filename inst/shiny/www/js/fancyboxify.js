$(document).ready( function() {
  
  // globals that will be used to save the width, height of the item that gets dblclicked
  var elem_width;
  var elem_height;
  var elem_parent;
  
  // the 'fancyboxify' function
  var fancyboxify = function() {
    
    // $this refers to the .gridster li element double-clicked
    $this = $(this).parent();
    
    // blur all other elements
    $("body > *").not($this).toggleClass("blur");
    
    // get the current browser dimensions
    var w = $(window).width();
    var h = $(window).height();
    var scrollTop = $(window).scrollTop();
    var scrollLeft = $(window).scrollLeft();
    
    var width = Math.min(1600, w-100);
    var height = Math.min(1200, h-100);
    
    var left = (w - width)/3 + scrollLeft;
    var top = (h - height)/3 + scrollTop;
    
    // get the height offset
    var height_offset = 0;
    $this.find("img").parent().siblings().each( function() {
      if (!$(this).hasClass("shiny-plot-output")) {
        // in reality we need to check borders, margins, padding,
        // but let's just add an extra 5px for each element we encounter...
        height_offset += parseFloat( $(this).css("height") ) + 5;
      }
    });
    
    // console.log("height_offset is: " + height_offset);
    
    // set the class of the double-clicked element
    if ($this.hasClass("fancybox")) {
      
      // display the constrols
      $("#gridster-control-container").css("display", "block");
      
      $this
        .find("div.shiny-plot-output")
        .css("width", elem_width)
        .css("height", elem_height)
      ;
      
      $this
        .toggleClass("fancybox")
        .removeAttr("style")
        .appendTo(elem_parent)
      ;
      
    } else {
      
      // hide the controls
      $("#gridster-control-container").css("display", "none");
      
      // update the 'reset' settings
      elem_width = $this
        .find("img")
        .attr("width")
      ;
      
      elem_height = $this
        .find("img")
        .attr("height")
      ;
      
      elem_parent = $this.parent();
      
      // fancyboxify
      $this
        .toggleClass("fancybox")
        .css("left", left + "px")
        .css("top", top + "px")
        .css("width", width)
        .css("height", height)
        .css("position", "absolute")
        .css("float", "left")
        .css("display", "inline")
        .css("padding", "20px")
        .css("background-color", $this.css("background-color"))
        .appendTo("body")
      ;
      
      $this
        .find("div.shiny-plot-output")
        .css("width", width)
        .css("height", height-height_offset)
      ;
      
      // disallow text highlighting
      $this
        .css("-webkit-user-select", "none")
        .css("-moz-user-select", "none")
        .css("-ms-user-select", "none")
      ;
      
    }
    
    // re-draw the plot
    $this
      .trigger("shown")
    ;
  };
  
  // double click to 'fancybox'ify
  // $(".gridster li").dblclick( fancyboxify );
  
  // click on a zoom icon to 'fancybox'ify
  $("i.icon-zoom-in").click( fancyboxify );
  
});
