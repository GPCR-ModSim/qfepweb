function destroy() {
  if (network !== null) {
    network.destroy();
    network = null;
  }
}

////////////////////////////////////////
// VISJS functions
////////////////////////////////////////
function draw(data) {
  destroy();
  nodes = [];
  edges = [];
  // create a network
  var container = document.getElementById("mynetwork");

  var options = {
  // ------------------------------------
  // GENERAL
   "nodes": {
   // could also make write PNGs as circular images, then set shape : image and omit size: n.
   // probably should do this. node clicking is inconsistent currently.
   //shape: "image",
    borderWidth: 1,
    shape: "image",
    shapeProperties: {
     useBorderWithImage : true
    },
    size: 120,
    font : {
     size : 40,
     color : 'blue'
    },
    color: {
      border: "#222222",
      background: "white",
      highlight: {
       border: '#2B7CE9',
       background: "white"
      }
    }
   },
   "edges": {
    "smooth": false,
    "font": {"size": 30},
   },
   // ------------------------------------
   // PHYSICS
   "physics": {
    "hierarchicalRepulsion": {
      "centralGravity": 0,
      "springLength": 0.3,
      "nodeDistance": 450,
      "damping": 0.13
    },
    "minVelocity": 0.76,
    "solver": "hierarchicalRepulsion",
    "timestep": 0.7
   },
   // ------------------------------------
   // MANIPULATION
   manipulation: {
    addNode : false,
    editEdge : false,
    addEdge: function (data, callback) {
     console.log('add edge', data);
     // empty if statement removes add edge to self funtionality.
     if (data.from == data.to) {
      // var r = confirm("Do you want to connect the node to itself?");
      // if (r === true) {
      //  callback(data);
      // }
     } else {
      callback(data);
     }
    }
  },
  interaction: {
   hover: true,
   multiselect: true
  }
 };

 network = new vis.Network(container, data, options);
 // disable physics after loading: https://stackoverflow.com/questions/32403578/stop-vis-js-physics-after-nodes-load-but-allow-drag-able-nodes
 network.on("stabilizationIterationsDone", function () {
  network.setOptions( { physics: false } );
 });

 // network.on("selectNode", function(ev){
 //  $("#network-popup").modal("show");
 //  $("h5.modal-title").text(ev.nodes[0]);
 // });
}

// function clearPopUp() {
//   document.getElementById("saveButton").onclick = null;
//   document.getElementById("cancelButton").onclick = null;
// }

// function cancelEdit(callback) {
//   clearPopUp();
//   callback(null);
// }

// function saveData(data, callback) {
//   data.id = document.getElementById("node-id").value;
//   data.label = document.getElementById("node-label").value;
//   clearPopUp();
//   callback(data);
// }
// SAVE SELECTION BUTTON
// saves a group of node IDs to a variable. Can be used downstream to look at statistics for the selected group.
function bs_alert(message, level){
 $('#alert-message').html(message);
 $('.alert').attr("class", "alert alert-dismissible " + level);
 $('.alert').show();
}

function saveselection(){

  const selNodes = data.nodes.get(network.getSelectedNodes());
  let labels = selNodes.map(n => n.label);

  if (selNodes.length == 0){
    // catch cases where user has not selected nodes.
    bs_alert("Please select ligands to save to group.", "alert-warning");
  } else {
  // require user to insert group name. Repop when no text input. :
  do{
    // prompt with string field
    var groupName = prompt("Enter an identifier for this group, then press OK. Press Cancel\n\
if you do not want to make a group with molecule(s):\n" + labels + ".");
} while(groupName !== null && groupName === "") // the while statement makes the prompt reappear
// when no string is input.

// if user has provided a name for the group, save the group name and node IDs.
if (groupName !== null) {
  console.log(groupName, selNodes)
  // HERE PROCESS groupName AND selNodes TO NEXT STEP IN WORKFLOW
  // TODO
  }}

}

// FINISH BUTTON/
// reduced version of https://stackoverflow.com/questions/40489700/visjs-save-manipulated-data-to-json

function savegraph(){
  // make variable that contains edges information.
  const nodePositions = data.edges.map(({ from, to }) => ({ from, to }))

  // save the edges as a JSON file for Q.
  writedata = JSON.stringify(nodePositions)

  bs_alert("Current network edges have been saved for Q." +
    "<br /> You can find a list of edges here (hlink).",
    "alert-info");

  // HERE PROCESS writedata TO NEXT STEP IN WORKFLOW
  // TODO
};



// Get out JSON data.
var network = null;

// FIXME: Take all nodes
function load_charge(charge=0){
  var nodes = new vis.DataSet(json_data[charge].nodes);
  var edges = new vis.DataSet(json_data[charge].edges);

  // create a network
  return {
    nodes: nodes,
    edges: edges
  };
}

////////////////////////////////////////
// undo/redo functions
////////////////////////////////////////

//initialize
// FIXME
//function history(){
//  let history_list_back = [];
//  let history_list_forward = [];
//
//  // initial data
//  history_list_back.push({
//    nodes_his: data.nodes.get(data.nodes.getIds()),
//    edges_his: data.edges.get(data.edges.getIds())
//  });
//  // event on
//  data.nodes.on("*", change_history_back);
//  data.edges.on("*", change_history_back);
//}
//
//function change_history_back() {
//  history_list_back.unshift({
//    nodes_his: data.nodes.get(data.nodes.getIds()),
//    edges_his: data.edges.get(data.edges.getIds())
//  });
//  //reset forward history
//  history_list_forward = [];
//  // apply css
//  css_for_manipulation();
//}
//function redo_css_active() {
//  $("#btn-undo").css({
//    "background-color": "inherit",
//    color: "#878787",
//    cursor: "pointer"
//  });
//};
//function undo_css_active() {
//  $("#btn-redo").css({
//    "background-color": "inherit",
//    color: "#878787",
//    cursor: "pointer"
//  });
//};
//
//function redo_css_inactive() {
//  $("#btn-undo").css({
//    "background-color": "inherit",
//    color: "#EBEBEB",
//    cursor: "inherit"
//  });
//};
//
//function undo_css_inactive() {
//  $("#btn-redo").css({
//    "background-color": "inherit",
//    color: "#EBEBEB",
//    cursor: "inherit"
//  });
//};
//
//function css_for_manipulation() {
//  if (history_list_back.length === 1) {
//    redo_css_inactive();
//  } else {
//    redo_css_active();
//  };
//  if (history_list_forward.length === 0) {
//    undo_css_inactive();
//  } else {
//    undo_css_active();
//  };
//};
//
//function undo(){
//  if (history_list_back.length > 1) {
//    const current_nodes = data.nodes.get(data.nodes.getIds());
//    const current_edges = data.edges.get(data.edges.getIds());
//    const previous_nodes = history_list_back[1].nodes_his;
//    const previous_edges = history_list_back[1].edges_his;
//    // event off
//    data.nodes.off("*", change_history_back);
//    data.edges.off("*", change_history_back);
//    // undo without events
//    if (current_nodes.length > previous_nodes.length) {
//      const previous_nodes_diff = _.differenceBy(
//        current_nodes,
//        previous_nodes,
//        "id"
//      );
//      data.nodes.remove(previous_nodes_diff);
//    } else {
//      data.nodes.update(previous_nodes);
//    }
//
//    if (current_edges.length > previous_edges.length) {
//      const previous_edges_diff = _.differenceBy(
//        current_edges,
//        previous_edges,
//        "id"
//      );
//      data.edges.remove(previous_edges_diff);
//    } else {
//      data.edges.update(previous_edges);
//    }
//    // recover event on
//    data.nodes.on("*", change_history_back);
//    data.edges.on("*", change_history_back);
//
//    history_list_forward.unshift({
//      nodes_his: history_list_back[0].nodes_his,
//      edges_his: history_list_back[0].edges_his
//    });
//    history_list_back.shift();
//            // apply css
//    css_for_manipulation();
//  }
//};
//
//function redo(){
//  if (history_list_forward.length > 0) {
//    const current_nodes = data.nodes.get(data.nodes.getIds());
//    const current_edges = data.edges.get(data.edges.getIds());
//    const forward_nodes = history_list_forward[0].nodes_his;
//    const forward_edges = history_list_forward[0].edges_his;
//    // event off
//    data.nodes.off("*", change_history_back);
//    data.edges.off("*", change_history_back);
//    // redo without events
//    if (current_nodes.length > forward_nodes.length) {
//      const forward_nodes_diff = _.differenceBy(
//        current_nodes,
//        forward_nodes,
//        "id"
//      );
//      data.nodes.remove(forward_nodes_diff);
//    } else {
//      data.nodes.update(forward_nodes);
//    }
//    if (current_edges.length > forward_edges.length) {
//      const forward_edges_diff = _.differenceBy(
//        current_edges,
//        forward_edges,
//        "id"
//      );
//      data.edges.remove(forward_edges_diff);
//    } else {
//      data.edges.update(forward_edges);
//    }
//    // recover event on
//    data.nodes.on("*", change_history_back);
//    data.edges.on("*", change_history_back);
//    history_list_back.unshift({
//      nodes_his: history_list_forward[0].nodes_his,
//      edges_his: history_list_forward[0].edges_his
//    });
//    // history_list_forward
//    history_list_forward.shift();
//        // apply css
//    css_for_manipulation();
//  }
//};

function shake(evt){
  network.setOptions(
    {physics:
      {enabled: !network.physics.physicsEnabled,
       timestep: 0.5}});
  if (network.physics.physicsEnabled) {
    evt.target.innerText = "Physics Off";
  }else{
    evt.target.innerText = "Physics On";
  }
}

function change_group(evt){
  nodes = new vis.DataSet(json_data[evt.target.value].nodes);
  edges = new vis.DataSet(json_data[evt.target.value].edges);

  // create a network
  data = {
    nodes: nodes,
    edges: edges
  };
  draw(data);
}

$(document).ready(function() {
              // apply css
  //css_for_manipulation();
  $("#save-sele").on("click", saveselection);
  //$("#btn-undo").on("click", undo);
  //$("#btn-redo").on("click", redo);
  $("#extract-positions").on("click", savegraph);
  $("#shake-up").on("click", shake);
  document.getElementById("selected-charge").addEventListener("change", change_group);
  var data = load_charge(document.getElementById("selected-charge").value);
  draw(data);
});
