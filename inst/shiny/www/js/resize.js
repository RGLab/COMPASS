$( function() {
  
  $(window).resize( function() {
    
    var width = $(window).width();
    var height = $(window).height();
    $(".gridster ul").gridster({
      widget_base_dimensions: [width, height]
    });
    
  });
  
});