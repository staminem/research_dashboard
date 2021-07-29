/*
$('#submit').click(function(event){
    // Prevent redirection with AJAX for contact form
    var form = $('#form');
    var form_id = 'form';
    var url = form.prop('action');
    var type = form.prop('method');
    var formData = getContactFormData(form_id);

    // submit form via AJAX
    send_form(form, form_id, url, type, modular_ajax, formData);
});

function getFormData(form) {
    // creates a FormData object and adds chips text
    var formData = new FormData(document.getElementById(form));
//    for (var [key, value] of formData.entries()) { console.log('formData', key, value);}
    return formData
}

function modular_ajax(url, type, formData) {
    // Most simple modular AJAX building block
    $.ajax({
        url: url,
        type: type,
        data: formData,
        processData: false,
        contentType: false,
        beforeSend: function() {
            // show the preloader (progress bar)
            $('#form-response').html("<div class='progress'><div class='indeterminate'></div></div>");
        },
        complete: function () {
            // hide the preloader (progress bar)
            $('#form-response').html("");
        },
        success: function ( data ){
            if ( !$.trim( data.feedback )) { // response from Flask is empty
                toast_error_msg = "An empty response was returned.";
                toast_category = "danger";
            }
            else { // response from Flask contains elements
                toast_error_msg = data.feedback;
                toast_category = data.category;
            }
        },
        error: function(xhr) {console.log("error. see details below.");
            console.log(xhr.status + ": " + xhr.responseText);
            toast_error_msg = "An error occured";
            toast_category = "danger";
        },
    }).done(function() {
        M.toast({html: toast_error_msg, classes: 'bg-' +toast_category+ ' text-white'});
    });
};

*/

/*
function processForm() {
    var userEmail = document.getElementById("user_email");
    var userTool = document.getElementById("user_tool");
    var pmids = document.getElementById("pmids[]");
    var keywords = document.getElementById("keywords[]");
    var serverData = [
        {"user_email": userEmail.value,
        "user_tool": userTool.value}
        ];

    console.log(userEmail.value);
    console.log(userTool.value);
    console.log(pmids.value);
    console.log(keywords.value);
    console.log("Hello World!")

    $.ajax({
        type: "POST",
        url: "/update_table",
        data: JSON.stringify(serverData),
        contentType: "application/json",
        dataType: "json",
        success: function(results){
            console.log(results);
        }

    });
}
*/

/*function store() {
  var pmids = document.getElementById('pmids').value;
  sessionStorage.setItem('pmids', pmids;
  }
*/

/*function retrieve() {
    var storedPmids = JSON.parse(sessionStorage.getItem("pmids"));
    if (storedPmids !== null) {
        var options = Array.from(document.getElementById("pmids[]").options);
        options.forEach(function(option) {
            for (var i = 0; i < storedPmids.length; i++) {
                if (option.value == storedPmids[i]) {
                    option.selected = true;
                }
            }
        });
    }
}*/

function spinner() {
    document.getElementById("img").style.display = "block";
}

function aradarPlot() {
    var radarPlot = document.getElementById("radarPlot");
    if (radarPlot.style.width == "80%"){
        radarPlot.style.width = "0%";
        radarPlot.style.height = "0px";
    }else{
        radarPlot.style.width = "80%";
        radarPlot.style.height = "800px";
    }
}

function abarPlot() {
    var barPlot = document.getElementById("barPlot");
    if (barPlot.style.display == "block"){
        barPlot.style.display = "none";
    }else{
        barPlot.style.display = "block";
    }
}

var pictures = ["../static/images/chart1.jpeg",
                "../static/images/chart2.jpeg",
                "../static/images/chart3.jpeg",
                "../static/images/chart4.jpeg",
                "../static/images/chart5.jpeg",
                "../static/images/chart6.jpeg",
                "../static/images/chart7.jpeg",
                "../static/images/chart8.jpeg",
                "../static/images/chart9.jpeg",
                "../static/images/chart10.jpeg",
                "../static/images/chart11.jpeg",
                "../static/images/chart12.jpeg",
                "../static/images/chart13.jpeg",
                "../static/images/chart14.jpeg",
                "../static/images/chart15.jpeg",
                ]

function showPic() {
    var randomNum = Math.floor(Math.random() * pictures.length);
    document.getElementById("picture").src = pictures[randomNum];
    }

function showModal() {
    $('#loadingModal').modal('toggle')
}
