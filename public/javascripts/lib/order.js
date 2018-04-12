'use strict';

let drawText = false;
let renderInterval = 33; //milliseconds
let gravity = false;
let gravConst = 32.2 * 12 / 1000000 / 8.25; // ft/s^2 -> heads/ms^2 (1 head = 8.25 in)
let frameRate = 0;
let collisions = 0;
let collideRate = 0;
let pause = false;

//[3] ... [0] help text
function ToolTip(now, x, y) {
    this.t = now;
    this.x = x;
    this.y = y;
}

ToolTip.prototype.render = function(context) {
    let e = Date.now() - this.t;
    let v;
    if (e < 500) {
        v = 255;
    } else if (e < 1500) {
        v = 255 * (1 - (e - 500) / 1000);
    } else {
        return false;
    }
    v = Number(Math.round(v)).toString(16);
    if (v.length === 1) {
        v = '0' + v;
    }
    let color = '#' + v + v + v;
    let y = this.y;
    for (let x = this.x; x <= this.x + 90; x+=90) {
        context.beginPath();
        context.strokeStyle = color;
        context.lineTo(x + 5, y);
        context.lineTo(x + 35, y);
        context.quadraticCurveTo(x + 40, y, x + 40, y + 5);
        context.lineTo(x + 40, y + 35);
        context.quadraticCurveTo(x + 40, y + 40, x + 35, y + 40);
        context.lineTo(x + 5, y + 40);
        context.quadraticCurveTo(x, y + 40, x, y + 35);
        context.lineTo(x, y + 5);
        context.quadraticCurveTo(x, y, x + 5, y);
        context.stroke();
    }

    context.beginPath();
    context.strokeStyle = color;
    context.fillStyle = color;
    context.arc(this.x + 55, this.y + 20, 2, 0, Math.PI*2);
    context.arc(this.x + 65, this.y + 20, 2, 0, Math.PI*2);
    context.arc(this.x + 75, this.y + 20, 2, 0, Math.PI*2);
    context.fill();

    context.textBaseline = 'top';
    context.font = '12pt Helvetica';
    context.fillText('3', this.x + 15, this.y + 10);
    context.fillText('0', this.x + 105, this.y + 10);

    return true;
}

let toolTip = null;

/*******************************
/*  Scene2D
/******************************/
function Scene2D(canvas) {
    this.canvas = canvas;
    this.matrix = new Matrix2D();
    this.matrixGood = false;
    this.context = canvas.getContext('2d');
    this.sceneWidth = canvas.width;
    this.sceneHeight = canvas.height;
    this.clipLeft = -this.sceneWidth*0.5;
    this.clipRight = this.sceneWidth*0.5;
    this.clipBottom = -this.sceneHeight*0.5;
    this.clipTop = this.sceneHeight*0.5;
    this.panX = 0;
    this.panY = 0;
    this.scale = 1.0;
    //this.rotation = 0;
    this.points = [];
    this.pointsTransformed = [];
    this.pointsLast = [];
    this.pointsAccel = [];
    this.pointsInvMass = [];
    this.numPoints = 0;
    this.polyItems = [];
    this.numPolyItems = 0;
    this.constraints = [];
    this.numConstraints = 0;
    this.dragItem = null;
    this.deleteDrag = false;
}

Scene2D.prototype.addPoint = function(x, y) {
    this.points[this.points.length] = this.pointsLast[this.pointsLast.length] = x;
    this.points[this.points.length] = this.pointsLast[this.pointsLast.length] = y;
    this.pointsAccel[this.pointsAccel.length] = 0;
    this.pointsAccel[this.pointsAccel.length] = 0;
    this.pointsInvMass[this.pointsInvMass.length] = 1;

    return this.numPoints++;
}

Scene2D.prototype.getPoint = function(idx) {
    let idx2D = idx*2;
    return [this.points[idx2D], this.points[idx2D+1]];
}

Scene2D.prototype.addPolyItem = function(item) {
    this.polyItems[this.polyItems.length] = item;

    return this.numPolyItems++;
}

Scene2D.prototype.makeSnapCurve = function(n, ll, ul) {
}

Scene2D.prototype.render = function() {
    let halfWidth = this.sceneWidth*0.5;
    let halfHeight = this.sceneHeight*0.5;

    if (!this.matrixGood) {
        this.matrix.identity();
        this.matrix.translate(halfWidth-this.panX, halfHeight-this.panY);
        this.matrix.scale(this.scale, this.scale);
        //this.matrix.rotate(this.rotation);
        this.matrixGood = true;
    }

    this.pointsTransformed = this.matrix.transformArray(this.points);

    this.context.save();
    this.context.fillStyle = '#000000';
    this.context.fillRect(0, 0, this.sceneWidth, this.sceneHeight);
    //this.context.fillStyle = '#333333';
    //this.context.fillRect(0, this.sceneHeight - 30, this.sceneWidth, 30);
    //this.context.fillStyle = '#aaaaaa';
    //this.context.textBaseline = 'top';
    //this.context.font = '12pt Helvetica';
    //this.context.fillText('joshua.l.salisbury@gmail.com', this.sceneWidth - 250, this.sceneHeight - 25);
    if (drawText) {
        this.context.fillStyle = '#aaaaaa';
        this.context.textBaseline = 'top';
        this.context.font = '12pt Helvetica';
        this.context.fillText('frame rate: ' + frameRate, 5, 5);
        this.context.fillText('collide rate: ' + collideRate + ' / sec', 5, 20);
        this.context.fillText('polys: ' + this.numPolyItems, 5, 35);
    }

    for (let i=0; i<this.numPolyItems; i++) {
        this.polyItems[i].render();
    }

    if (toolTip !== null) {
        if (toolTip.render(this.context) === false) {
            toolTip = null;
        }
    }

    //if (this.dragItem != null) {
    //    this.dragItem.render(this);
    //}

    this.context.restore();
}

Scene2D.prototype.accumulateForces = function() {
    for (let i=0; i<this.numPoints; i++) {
        let i2 = i*2;
        let xIdx = i2;
        let yIdx = i2+1;

        //this.points3DAccel[xIdx] = gravity?gravConst*windCos*0.1:0;
        //this.points3DAccel[yIdx] = gravity?gravConst*windSin*0.1:0;
        this.pointsAccel[xIdx] = 0;
        this.pointsAccel[yIdx] = gravity?-gravConst:0; // 1g of gravity in heads (8.25in: female) per millisec^2
    }

    this.snapPolys();
}

Scene2D.prototype.verlet = function(timeDel) {
    for (let i=0; i<this.numPoints; i++) {
        let i2 = i*2;
        let xIdx = i2;
        let yIdx = i2+1;
        let td2 = timeDel*timeDel;

        let tempX = this.points[xIdx];
        let tempY = this.points[yIdx];

        //this.points[xIdx] += (this.points[xIdx]-this.pointsLast[xIdx])*0.99 + this.pointsAccel[xIdx]*timeDel*timeDel;
        //this.points[yIdx] += (this.points[yIdx]-this.pointsLast[yIdx])*0.99 + this.pointsAccel[yIdx]*timeDel*timeDel;
        this.points[xIdx] += this.points[xIdx]-this.pointsLast[xIdx] + this.pointsAccel[xIdx]*td2;
        this.points[yIdx] += this.points[yIdx]-this.pointsLast[yIdx] + this.pointsAccel[yIdx]*td2;

        this.pointsLast[xIdx] = tempX;
        this.pointsLast[yIdx] = tempY;
    }
}

Scene2D.prototype.constrain = function() {
    for (let i=0; i<this.numPolyItems; i++) {
        this.polyItems[i].constrain();
    }
    if (this.dragItem != null) {
        this.dragItem.constrain();
    }
}

Scene2D.prototype.clip = function() {
    let temp = 0;
    for (let i=0; i<this.numPoints; i++) {
        if (this.points[i*2] < this.clipLeft) {
            this.points[i*2] = this.clipLeft;
        }
        if (this.clipRight < this.points[i*2]) {
            this.points[i*2] = this.clipRight;
        }
        if (this.points[i*2+1] < this.clipBottom) {
            this.points[i*2+1] = this.clipBottom;
        }
        if (this.clipTop < this.points[i*2+1]) {
            this.points[i*2+1] = this.clipTop;
        }
    }
}

Scene2D.prototype.collidePolys = function() {
    //for (let i=0; i<this.numPolyItems; i++) {
    //  this.polyItems[i].color = '#eeeeee';
    //}
    let collision;

    for (let i=0; i<this.numPolyItems; i++) {
        for (let j=i+1; j<this.numPolyItems; j++) {
            if (this.polyItems[i].interacting && this.polyItems[j].interacting) {
                collision = this.polyItems[i].collidePoly(this.polyItems[j]);
                if (collision) {
                    collisions++;
                }
                //this.polyItems[i].color = '#ee0000';
                //this.polyItems[j].color = '#ee0000';
            }
        }
    }
}

Scene2D.prototype.snapPolys = function() {
    //for (let i=0; i<this.numPolyItems; i++) {
    //  this.polyItems[i].color = '#eeeeee';
    //}
    for (let i=0; i<this.numPolyItems; i++) {
        for (let j=i+1; j<this.numPolyItems; j++) {
            this.polyItems[i].snapPoly(this.polyItems[j]);
            //if (this.polyItems[i].snapPoly(this.polyItems[j])) {
            //  this.polyItems[i].color = '#ee0000';
            //  this.polyItems[j].color = '#ee0000';
            //}
        }
    }
}

Scene2D.prototype.fixPolys = function() {
    for (let i=0; i<this.numPolyItems; i++) {
        if (this.polyItems[i].isInverted()) {
            this.polyItems[i].invert();
        }
    }
}

Scene2D.prototype.safeActions = function() {
    if (this.deleteDrag) {
        if (this.dragItem) {
            this.polyItems[this.dragItem.polyItemIdx].interacting = true;
            this.dragItem = null;
        }
        this.deleteDrag = false;
    }
    this.fixPolys();
}

Scene2D.prototype.tic = function() {
    this.render();
    this.accumulateForces();
    this.verlet(renderInterval);
    this.constrain();
    this.collidePolys();
    this.collidePolys();
    this.constrain();
    this.collidePolys();
    this.collidePolys();
    this.clip();
    this.safeActions();
}

Scene2D.prototype.resize = function() {
    this.sceneWidth = this.canvas.width;
    this.sceneHeight = this.canvas.height;
    this.clipLeft = -this.sceneWidth*0.5;
    this.clipRight = this.sceneWidth*0.5;
    this.clipBottom = -this.sceneHeight*0.5;
    this.clipTop = this.sceneHeight*0.5;
    this.matrixGood = false;
}
/*******************************
/*  END Scene2D
/******************************/

/*******************************
/*  PolyItem
/******************************/
function PolyItem(scene, pIdx) {
    this.scene = scene;
    this.pIdx = pIdx;
    this.p2Didx = [];
    this.interacting = true;

    let i, j, n, v;

    for (i=0; i<this.pIdx.length; i++) {
        this.p2Didx[i] = 2 * this.pIdx[i];
    }
    this.sideLengths = [];
    this.invSideLengths = [];
    this.restLengths = [];
    this.restLengthsSq = [];
    this.edgeVecIdx = [];
    this.normalDepths = [];
    this.color = '#eeeeee';

    // INIT
    let points = this.scene.points;
    let minSelfProj = 0;
    let norm, proj;
    for (i=0; i<this.p2Didx.length; i++) {
        let si = (i+1<this.p2Didx.length)?i+1:0; // successor of i

        // compute side lengths
        v = [points[this.p2Didx[si]]-points[this.p2Didx[i]], points[this.p2Didx[si]+1]-points[this.p2Didx[i]+1]];
        this.sideLengths[i] = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
        this.invSideLengths[i] = 1.0 / this.sideLengths[i];

        // constraint length pairs
        this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[i];
        this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[si];
        n = this.p2Didx.length-3-(i>1?i-1:0);
        for (let j=i+2; j<i+n+2; j++) {
            this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[i];
            this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[j];
        }

        // depth from each face
        norm = [(points[this.p2Didx[si]+1]-points[this.p2Didx[i]+1])*this.invSideLengths[i], // right-hand unit normal of i-th side
                    (points[this.p2Didx[i]]-points[this.p2Didx[si]])*this.invSideLengths[i]];
        this.normalDepths[i] = 0;
        for (j=0; j<this.p2Didx.length; j++) {
            proj = (points[this.p2Didx[j]]-points[this.p2Didx[i]])*norm[0] + (points[this.p2Didx[j]+1]-points[this.p2Didx[i]+1])*norm[1];
            this.normalDepths[i] = (this.normalDepths[i] < -proj)?-proj:this.normalDepths[i];
        }
    }

    // constraint rest lengths
    for (i=0; i<this.edgeVecIdx.length/2; i++) {
        v = [points[this.edgeVecIdx[i*2+1]]-points[this.edgeVecIdx[i*2]], points[this.edgeVecIdx[i*2+1]+1]-points[this.edgeVecIdx[i*2]+1]];
        this.restLengths[i] = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
        this.restLengthsSq[i] = this.restLengths[i]*this.restLengths[i];
    }
}

PolyItem.prototype.render = function() {
    let points = this.scene.pointsTransformed;
    let context = this.scene.context;

    context.beginPath();
    context.strokeStyle = this.color;
    context.lineWidth = 3;
    //context.fillStyle = '#eeeeee';
    for (let i=0; i<this.p2Didx.length; i++) {
        context.lineTo(points[this.p2Didx[i]], this.scene.sceneHeight - points[this.p2Didx[i]+1]);
        //context.fillText(i, points[this.p2Didx[i]], this.scene.sceneHeight - points[this.p2Didx[i]+1]);
    }
    context.lineTo(points[this.p2Didx[0]], this.scene.sceneHeight - points[this.p2Didx[0]+1]); // close the polygon
    context.stroke();
}

PolyItem.prototype.constrain = function() {
    let points = this.scene.points;
    let vecH, delta, del, delSq, im1, im2, diff;
    let i, j;
    for (i=0; i<this.edgeVecIdx.length/2; i++) {
        vecH = i*2;
        delta = [points[this.edgeVecIdx[vecH+1]]-points[this.edgeVecIdx[vecH]],
                     points[this.edgeVecIdx[vecH+1]+1]-points[this.edgeVecIdx[vecH]+1]];
        delSq = delta[0]*delta[0] + delta[1]*delta[1];
        im1 = 1;//this.scene.pointsInvMass[this.edgeVecIdx[vecH]/2];
        im2 = 1;//this.scene.pointsInvMass[this.edgeVecIdx[vecH+1]/2];
        diff = (delSq-this.restLengthsSq[i])/((delSq+this.restLengthsSq[i])*(im1+im2)*2);
        for (j=0; j<2; j++) {
            del = diff*delta[j];
            points[this.edgeVecIdx[vecH]+j]   += del*im1;
            points[this.edgeVecIdx[vecH+1]+j] -= del*im2;
        }
    }
}

PolyItem.prototype.collidePoly = function(poly2) {
    let points = this.scene.points;

    let minOverlap = 1000;
    let bestNorm = [];
    let besti = -1;
    let bestj = -1;
    let bestp = -1;

    let i, j, p;
    let norm, proj;
    let min12Proj, max12Proj, min12Projp;
    let min21Proj, max21Proj, min21Projp;
    let minProj, maxProj, overlap

    let ijHitF, l, iHitP, jHitP, push;

    // First check this poly's sides for penetration against "that" poly's points
    for (i=0; i<this.p2Didx.length; i++) {
        j = (i+1<this.p2Didx.length)?i+1:0;
        norm = [(points[this.p2Didx[j]+1]-points[this.p2Didx[i]+1])*this.invSideLengths[i], // unit normal of i-th side
                    (points[this.p2Didx[i]]-points[this.p2Didx[j]])*this.invSideLengths[i]];

        min12Proj = 0;
        max12Proj = 0;
        min12Projp = -1;
        for (p=0; p<poly2.p2Didx.length; p++) {
            proj = (points[poly2.p2Didx[p]]-points[this.p2Didx[i]])*norm[0] + (points[poly2.p2Didx[p]+1]-points[this.p2Didx[i]+1])*norm[1];
            if (proj<min12Proj) {
                min12Proj = proj;
                min12Projp = p;
            }
            max12Proj = (max12Proj<proj)?proj:max12Proj;
        }

        if (min12Projp == -1) {
            return false;
        }

        minProj = (-this.normalDepths[i]<min12Proj)?-this.normalDepths[i]:min12Proj;
        maxProj = (max12Proj>0)?max12Proj:0;

        overlap = max12Proj - min12Proj + this.normalDepths[i] - maxProj + minProj;
        if (overlap > 0) {
            if (overlap < minOverlap) {
                minOverlap = overlap;
                bestNorm = norm;
                besti = this.p2Didx[i];
                bestj = this.p2Didx[j];
                bestp = poly2.p2Didx[min12Projp];
            }
        } else {
            return false;
        }
    }

    // Then check "that" poly's sides
    for (i=0; i<poly2.p2Didx.length; i++) {
        j = (i+1<poly2.p2Didx.length)?i+1:0;
        norm = [(points[poly2.p2Didx[j]+1]-points[poly2.p2Didx[i]+1])*poly2.invSideLengths[i], // unit normal of i-th side
                    (points[poly2.p2Didx[i]]-points[poly2.p2Didx[j]])*poly2.invSideLengths[i]];

        min21Proj = 0;
        max21Proj = 0;
        min21Projp = -1;
        for (p=0; p<this.p2Didx.length; p++) {
            proj = (points[this.p2Didx[p]]-points[poly2.p2Didx[i]])*norm[0] + (points[this.p2Didx[p]+1]-points[poly2.p2Didx[i]+1])*norm[1];
            if (proj<min21Proj) {
                min21Proj = proj;
                min21Projp = p;
            }
            max21Proj = (max21Proj<proj)?proj:max21Proj;
        }

        if (min21Projp == -1) {
            return false;
        }

        minProj = (-poly2.normalDepths[i]<min21Proj)?-poly2.normalDepths[i]:min21Proj;
        maxProj = (max21Proj>0)?max21Proj:0;

        overlap = max21Proj - min21Proj + poly2.normalDepths[i] - maxProj + minProj;
        if (overlap > 0) {
            //alert(overlap);
            if (overlap < minOverlap) {
                minOverlap = overlap;
                bestNorm = norm;
                besti = poly2.p2Didx[i];
                bestj = poly2.p2Didx[j];
                bestp = this.p2Didx[min21Projp];
            }
        } else {
            return false;
        }
    }

    // Update points
    if (Math.abs(points[bestj] - points[besti]) > Math.abs(points[bestj+1] - points[besti+1])) {
        ijHitF = Math.abs(points[bestp]-points[besti]) / Math.abs(points[bestj]-points[besti]);
    } else {
        ijHitF = Math.abs(points[bestp+1]-points[besti+1]) / Math.abs(points[bestj+1]-points[besti+1]);
    }
    l = 1.0/(ijHitF*ijHitF + (1-ijHitF)*(1-ijHitF));
    iHitP = (1-ijHitF)*l;
    jHitP = ijHitF*l;
    push = [bestNorm[0]*minOverlap*0.5, bestNorm[1]*minOverlap*0.5];
    //alert (bestNorm[0] + ' ' + bestNorm[1]);
    points[besti]   -= push[0]*iHitP;
    points[besti+1] -= push[1]*iHitP;
    points[bestj]   -= push[0]*jHitP;
    points[bestj+1] -= push[1]*jHitP;
    points[bestp]   += push[0];
    points[bestp+1] += push[1];

    return true;
}

PolyItem.prototype.snapPoly = function(poly2) {
    var points = this.scene.points;
    var pointsLast = this.scene.pointsLast;
    var pointsAccel = this.scene.pointsAccel;

    // check all this poly's sides against all that poly's sides
    for (var i=0; i<this.p2Didx.length; i++) {
        var si = (i+1<this.p2Didx.length)?i+1:0;

        for (var j=0; j<poly2.p2Didx.length; j++) {
            var sj = (j+1<poly2.p2Didx.length)?j+1:0;

            var v1 = [points[poly2.p2Didx[sj]]-points[this.p2Didx[i]], points[poly2.p2Didx[sj]+1]-points[this.p2Didx[i]+1]];
            var v2 = [points[poly2.p2Didx[j]]-points[this.p2Didx[si]], points[poly2.p2Didx[j]+1]-points[this.p2Didx[si]+1]];
            var ds1 = v1[0]*v1[0]+v1[1]*v1[1]; // square of distance from this poly's point i to that poly's point sj
            var ds2 = v2[0]*v2[0]+v2[1]*v2[1];

            if (ds1<100.0 && ds2<100.0) {
                var snapGap = 0.01;
                var ui = [(points[this.p2Didx[si]]-points[this.p2Didx[i]])*this.invSideLengths[i], (points[this.p2Didx[si]+1]-points[this.p2Didx[i]+1])*this.invSideLengths[i]];
                var uj = [(points[poly2.p2Didx[sj]]-points[poly2.p2Didx[j]])*poly2.invSideLengths[j], (points[poly2.p2Didx[sj]+1]-points[poly2.p2Didx[j]+1])*poly2.invSideLengths[j]];
                var ni = [ui[1], -ui[0]];
                var nj = [uj[1], -uj[0]];
                var v1 = [points[poly2.p2Didx[sj]]+nj[0]*snapGap-points[this.p2Didx[i]], points[poly2.p2Didx[sj]+1]+nj[1]*snapGap-points[this.p2Didx[i]+1]];
                var v2 = [points[poly2.p2Didx[j]]-points[this.p2Didx[si]]+nj[0]*snapGap, points[poly2.p2Didx[j]+1]-points[this.p2Didx[si]+1]+nj[1]*snapGap];

                var vi = [points[this.p2Didx[i]]-pointsLast[this.p2Didx[i]], points[this.p2Didx[i]+1]-pointsLast[this.p2Didx[i]+1]];
                var vsi = [points[this.p2Didx[si]]-pointsLast[this.p2Didx[si]], points[this.p2Didx[si]+1]-pointsLast[this.p2Didx[si]+1]];
                var vj = [points[poly2.p2Didx[j]]-pointsLast[poly2.p2Didx[j]], points[poly2.p2Didx[j]+1]-pointsLast[poly2.p2Didx[j]+1]];
                var vsj = [points[poly2.p2Didx[sj]]-pointsLast[poly2.p2Didx[sj]], points[poly2.p2Didx[sj]+1]-pointsLast[poly2.p2Didx[sj]+1]];

                pointsAccel[this.p2Didx[i]]     +=  0.0001*v1[0] - vi[0]*0.0001;
                pointsAccel[this.p2Didx[i]+1]   +=  0.0001*v1[1] - vi[1]*0.0001;
                pointsAccel[this.p2Didx[si]]    +=  0.0001*v2[0] - vsi[0]*0.0001;
                pointsAccel[this.p2Didx[si]+1]  +=  0.0001*v2[1] - vsi[1]*0.0001;
                pointsAccel[poly2.p2Didx[j]]    += -0.0001*v2[0] - vj[0]*0.0001;
                pointsAccel[poly2.p2Didx[j]+1]  += -0.0001*v2[1] - vj[1]*0.0001;
                pointsAccel[poly2.p2Didx[sj]]   += -0.0001*v1[0] - vsj[0]*0.0001;
                pointsAccel[poly2.p2Didx[sj]+1] += -0.0001*v1[1] - vsj[1]*0.0001;

                return true;
            }
        }
    }

    return false;
}

PolyItem.prototype.pointIntersects = function(point) {
    let points = this.scene.points;
    let i, j, norm, proj;

    for (i=0; i<this.p2Didx.length; i++) {
        j = (i+1<this.p2Didx.length)?i+1:0;
        norm = [(points[this.p2Didx[j]+1]-points[this.p2Didx[i]+1])*this.invSideLengths[i], // unit normal of i-th side
                    (points[this.p2Didx[i]]-points[this.p2Didx[j]])*this.invSideLengths[i]];
        proj = (point[0]-points[this.p2Didx[i]])*norm[0] + (point[1]-points[this.p2Didx[i]+1])*norm[1];
        if (proj > 0) {
            return false;
        }
    }

    return true;
}

// checks the cross product of the first two sides. this assumes a convex poly
PolyItem.prototype.isInverted = function() {
    let points = this.scene.points;
    let v1 = [points[this.p2Didx[1]] - points[this.p2Didx[0]], points[this.p2Didx[1]+1] - points[this.p2Didx[0]+1]];
    let v2 = [points[this.p2Didx[2]] - points[this.p2Didx[1]], points[this.p2Didx[2]+1] - points[this.p2Didx[1]+1]];
    return v1[0]*v2[1]-v1[1]*v2[0] < 0;
}

// reverses the vertices' direction
PolyItem.prototype.invert = function() {
    let inverted = [];
    for (let i=this.p2Didx.length-1; i>=0; i--) {
        inverted[inverted.length] = this.p2Didx[i];
    }
    this.p2Didx = inverted;
}
/*******************************
/*  END PolyItem
/******************************/

/*******************************
/*  DragItem
/******************************/
function DragItem(scene, dragP, polyItemIdx) {
    this.scene = scene;
    this.dragP = dragP;
    this.polyItemIdx = polyItemIdx;
    this.pIdx = scene.polyItems[polyItemIdx].pIdx;
    this.p2Didx = [];
    let i, v;

    for (i=0; i<this.pIdx.length; i++) {
        this.p2Didx[i] = 2 * this.pIdx[i];
    }
    this.restLengths = [];
    this.restLengthsSq = [];
    this.color = '#eeeeee';

    // INIT
    let points = this.scene.points;
    // constraint rest lengths
    for (i=0; i<this.p2Didx.length; i++) {
        v = [points[this.p2Didx[i]]-this.dragP[0], points[this.p2Didx[i]+1]-this.dragP[1]];
        this.restLengths[i] = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
        this.restLengthsSq[i] = this.restLengths[i]*this.restLengths[i];
    }
}

DragItem.prototype.render = function() {
    let points = this.scene.pointsTransformed;
    let context = this.scene.context;
    let dragPointTransformed = this.scene.matrix.transformArray(this.dragP);

    context.beginPath();
    context.strokeStyle = this.color;
    for (let i=0; i<this.p2Didx.length; i++) {
        context.moveTo(dragPointTransformed[0], this.scene.sceneHeight - dragPointTransformed[1]);
        context.lineTo(points[this.p2Didx[i]], this.scene.sceneHeight - points[this.p2Didx[i]+1]);
    }
    context.stroke();
}

DragItem.prototype.constrain = function() {
    let points = this.scene.points;
    let i, j;
    let delta, delSq, diff, del;

    for (i=0; i<this.p2Didx.length; i++) {
        delta = [points[this.p2Didx[i]]-this.dragP[0],
                     points[this.p2Didx[i]+1]-this.dragP[1]]
        delSq = delta[0]*delta[0] + delta[1]*delta[1];
        diff = (delSq-this.restLengthsSq[i])/(delSq+this.restLengthsSq[i]);
        for (j=0; j<2; j++) {
            del = diff*delta[j];
            points[this.p2Didx[i]+j] -= del;
        }
    }
}
/*******************************
/*  END DragItem
/******************************/

/*******************************
/*  MAIN
/******************************/
let baseCanvas = document.getElementById('baseCanvas');
let scene = new Scene2D(baseCanvas);
let state = 'none';
let px = 0;
let py = 0;
let i;

let colorMap = {
    10: '#ff0000', // red
    9:  '#ff9d00', // orange
    8:  '#fffb00', // yellow
    7:  '#80ff00', // green
    6:  '#00c8ff', // cyan
    5:  '#0040ff', // blue
    4:  '#8000ff', // indigo
    3:  '#dd00ff'  // fuscia
}

function genPoly(n, x, y) {
    let t = 2*Math.PI/n
    let sl = Math.sin(t/2);
    let r = 50 / sl;

    let p = [];
    for (let i=0; i<n; i++) {
        p[2*i] = x + r * Math.cos(t/2 + i*t);
        p[2*i+1] = y + r * Math.sin(t/2 + i*t);
    }

    return p;
}

function addPoly(n) {
    //let x = Math.random() * baseCanvas.width() - baseCanvas.width()/2;
    //let y = Math.random() * baseCanvas.height() - baseCanvas.height()/2;
    let p = genPoly(n, px, py);
    let pIdx = [];
    for (let i = 0; i < p.length; i+=2) {
        pIdx.push(scene.addPoint(p[i], p[i+1]));
    }
    let polyItem = new PolyItem(scene, pIdx);
    polyItem.color = colorMap[n];
    scene.addPolyItem(polyItem);
}

function addRhombA() {
    let l = 50.0;
    let p1 = scene.addPoint(px+l*Math.cos(Math.PI/10), py);
    let p2 = scene.addPoint(px, py+l*Math.sin(Math.PI/10));
    let p3 = scene.addPoint(px-l*Math.cos(Math.PI/10), py);
    let p4 = scene.addPoint(px, py-l*Math.sin(Math.PI/10));

    scene.addPolyItem(new PolyItem(scene, [p1,p2,p3,p4]));
}

function addRhombB() {
    let l = 50.0;
    let p1 = scene.addPoint(px+l*Math.cos(Math.PI/5), py);
    let p2 = scene.addPoint(px, py+l*Math.sin(Math.PI/5));
    let p3 = scene.addPoint(px-l*Math.cos(Math.PI/5), py);
    let p4 = scene.addPoint(px, py-l*Math.sin(Math.PI/5));

    scene.addPolyItem(new PolyItem(scene, [p1,p2,p3,p4]));
}

let frames=0;

setInterval(function () {
    scene.tic();
    frames++;
}, renderInterval);

setInterval(function () {
    frameRate = frames;
    collideRate = collisions;
    frames=0;
    collisions=0;
}, 1000);

document.addEventListener("keydown", function(event) {
    if (state == 'alert_next_key') {
        state = 'none';
        alert(event.which);
    }

    // KEYS
    switch(event.which) {
        case 8:                                       //backspace
          break;
        case 9:                                       //tab
          break;
        case 13:                                      //enter
          break;
        case 27:                                      //esc
          break;
        case 37:                                      //left arrow
          if (event.shiftKey) {
          } else {
              scene.panX -= 5;
          }
          break;
        case 38:                                      //up arrow
          scene.panY -= 5;
          break;
        case 39:                                      //right arrow
          if (event.shiftKey) {
          } else {
              scene.panX += 5;
          }
          break;
        case 40:                                      //down arrow
          scene.panY += 5;
          break;
        case 46:                                      //delete
          break;
        case 48:                                      //'0'
          addPoly(10)
          break;
        case 49: case 50:                             //'1'-'2'
          break;
        case 51: case 52: case 53: case 54: case 55: case 56: case 57: //'3'-'9'
          addPoly(event.which - 48);
          break;
        case 61:                                      //'=' / '+'
          break;
        case 65:                                      //'a'
          addRhombA();
          break;
        case 66:                                      //'b'
          addRhombB();
          break;
        case 67:                                      //'c'
          break;
        case 68:                                      //'d'
          addDart(); //TODO
          break;
        case 69:                                      //'e'
          break;
        case 70:                                      //'f'
          addPoly(Math.floor(Math.random()*6+3));
          break;
        case 71:                                      //'g'
          gravity = !gravity;
          break;
        case 72:                                      //'h'
          break;
        case 73:                                      //'i'
          scene.scale *= 1.1;
          break;
        case 76:                                      //'l'
          break;
        case 78:                                      //'n'
          drawText = !drawText;
          break;
        case 79:                                      //'o'
          scene.scale /= 1.1;
          break;
        case 80:                                      //'p'
          pause = !pause;
          break;
        case 82:                                      //'r'
          break;
        case 83:                                      //'s'
          break;
        case 85:                                      //'u'
          scene.zTranslate -= 0.2;
          break;
        case 86:                                      //'v'
          break;
        case 87:                                      //'w'
          break;
        case 89:                                      //'y'
          break;
        case 90:                                      //'z'
          break;
        case 173:                                     //'-'
          break;
        case 190:                                     //'.'
          break;
        case 191:                                     //'/'
          break;
        case 192:                                     //'`'
          if (event.shiftKey) {
              state = 'alert_next_key';
          }
          break;
        case 220:                                     //'\'
          break;
    }
});

baseCanvas.addEventListener('mousedown', function(event) {
    if (scene.numPolyItems === 0) {
        let canvasBound = baseCanvas.getBoundingClientRect();
        toolTip = new ToolTip(Date.now(), event.clientX - canvasBound.left + 10, event.clientY + canvasBound.top + 15);
    }
    for (let i=0; i<scene.numPolyItems; i++) {
        if (scene.polyItems[i].pointIntersects([px, py])) {
            scene.polyItems[i].interacting = false;
            scene.dragItem = new DragItem(scene, [px, py], i);
            break;
        }
    }
});

baseCanvas.addEventListener('mousemove', function(event) {
    let canvasBound = baseCanvas.getBoundingClientRect();
    px = -canvasBound.width * 0.5 + event.clientX - canvasBound.left;
    py = canvasBound.height * 0.5 - event.clientY + canvasBound.top;
    if (scene.dragItem != null) {
        scene.dragItem.dragP = [px, py];
    }
});

baseCanvas.addEventListener('mouseup', function(event) {
    scene.deleteDrag = true;
});

window.onresize = function() {
    baseCanvas.width = window.innerWidth;
    baseCanvas.height = window.innerHeight;
    scene.resize();
}

let event = new Event('resize');
window.dispatchEvent(event);
