[

   {"selector":"node", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "background-color": "white",
       "border-color": "black",
       "shape": "ellipse",
       "width": 60,
       "height": 60,
       "content": "data(id)",
       "border-width": "1px"
       }},
       
    {"selector":"node[kld>0]", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "background-color": "mapData(kld, 0, 1, white, orangered)"
       }},
       
    {"selector":"node[kld<0]", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "background-color": "mapData(kld, -1, 0, dodgerblue, white)"
       }},

    {"selector":"node[type=0]:", "css": {
       "shape": "round-rectangle"
       }},

    {"selector":"node[type=1]:", "css": {
       "shape": "ellipse"
       }},
       
    {"selector":"node[absorbed='y']:", "css": {
       "border-color": "red"
       }},
       
    {"selector":"node[absorbed='n']:", "css": {
       "border-color": "black"
       }},

   {"selector":"node:selected", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "border-color": "black",
       "content": "data(id)",
       "border-width": "3px",
       "overlay-opacity": 0.5,
       "overlay-color": "blue"
       }},

    {"selector": "edge", "css": {
        "line-color": "maroon",
        "source-arrow-color": "black",
        "target-arrow-color": "black",
        "target-arrow-shape": "triangle",
        "curve-style": "bezier"
        }},

    {"selector": "edge[score<=0]", "css": {
        "line-color": "mapData(score, -30, 0, red, lightGray)",
        "source-arrow-shape": "circle",
        "source-arrow-color": "orange",
        "target-arrow-shape": "tee",
        "target-arrow-color": "black",
        "curve-style": "bezier"
        }},

    {"selector": "edge[score>0]", "css": {
        "line-color": "mapData(score, 0, 30, lightGray, green)",
        "source-arrow-shape": "circle",
        "source-arrow-color": "orange",
        "target-arrow-shape": "tee",
        "target-arrow-color": "black",
        "curve-style": "bezier"
        }},

    {"selector": "edge:selected", "css": {
       "overlay-opacity": 0.2,
       "overlay-color": "maroon"
        }}

   ]
