'use strict';

let canvas = document.getElementById('baseCanvas');
let context = canvas.getContext('2d');

let points = [];
let points_last = [];
let points_accel = [];
let num_points = 0;
let points_transformed = [];

let poly_items = [];
let drag_item = null;
let delete_drag = false;

let palette = {
    width: -3,
    handle_width: 50,
    off_screen_width: 100,
    grip: null,
    clipping_plane: make_clipping_plane(
        -canvas.width * 0.5,
        -canvas.width * 0.5 + this.width,
        -canvas.height * 0.5,
        canvas.height * 0.5
    ),

    resize(width) {
        width = width < -3 ? -3 : width;
        width = width > 300 ? 300 : width;
        this.width = width;
        this.off_screen_width = (100 - width < 0) ? 0 : 100 - width;

        this.clipping_plane.clip_left = -canvas.width * 0.5 - this.off_screen_width;
        this.clipping_plane.clip_right = -canvas.width * 0.5 + this.width;
        this.clipping_plane.clip_bottom = -canvas.height * 0.5;
        this.clipping_plane.clip_top = canvas.height * 0.5;
    },

    render() {
        let points = [
            -canvas.width * 0.5,                  -canvas.height * 0.5,     // corner of palette
            -canvas.width * 0.5 + this.width,     -canvas.height * 0.5,     // corner of palette handle
            -canvas.width * 0.5 + this.width + 8, -canvas.height * 0.5 + 5, // corner of '+' text
            -canvas.width * 0.5 + this.width + 8, -canvas.height * 0.5 + 65 // corner of '+' text
        ];
        let points_transformed = matrix.transformArray(points);

        context.fillStyle = '#' + rgba_to_hex(51, 51, 51, 128);
        context.fillRect(points_transformed[0], points_transformed[1], this.width, canvas.height);
        context.fillStyle = '#' + rgba_to_hex(102, 102, 102, 128);
        context.fillRect(points_transformed[2], points_transformed[3], this.handle_width, canvas.height);

        context.fillStyle = '#aaaaaa';
        context.textBaseline = 'top';
        context.font = '60px Helvetica';
        context.fillText('\u002b', points_transformed[4], points_transformed[5]); // 'plus'
        context.fillText('\u2212', points_transformed[6], points_transformed[7]); // 'minus'
    },

    does_handle_intersect(point) {
        return point[0] > (-canvas.width * 0.5 + this.width) && point[0] < (-canvas.width * 0.5 + this.width + this.handle_width);
    },

    get_handle_point(point) {
        return [
            point[0] - this.width,
            point[1]
        ];
    },

    does_plus_intersect(point) {
        return point[0] > (-canvas.width * 0.5 + this.width)
            && point[0] < (-canvas.width * 0.5 + this.width + this.handle_width)
            && point[1] > (-canvas.height * 0.5)
            && point[1] < (-canvas.height * 0.5 + 60); // 60px Helvetica
    },

    double_size() {
        // TODO add a function as a property to poly items that generates their coordinates
    },

    does_minus_intersect(point) {
        return point[0] > (-canvas.width * 0.5 + this.width)
            && point[0] < (-canvas.width * 0.5 + this.width + this.handle_width)
            && point[1] > (-canvas.height * 0.5 + 60) // 60px Helvetica
            && point[1] < (-canvas.height * 0.5 + 120);
    },

    does_body_intersect(point) {
        return point[0] < this.width;
    }
};

let matrix = new Matrix2D();
let matrix_good = false;

let main_clipping_plane = make_clipping_plane(
    -canvas.width * 0.5,
    canvas.width * 0.5,
    -canvas.height * 0.5,
    canvas.height * 0.5
);

let clipping_planes = [
    main_clipping_plane,
    palette.clipping_plane
];

let draw_text = false;
let render_interval = 33; //milliseconds
let pan_x = 0;
let pan_y = 0;
let scale = 1.0;
let gravity_constant = 32.2 * 12 / 1000000;

let frames = 0;
let frame_rate = 0;
let collisions = 0;
let collide_rate = 0;

let state = 'none';
let px = 0, py = 0;

function make_poly_item(point_generator, x, y, num_sides, side_length) {
    let p = point_generator(x, y, num_sides, side_length);

    let p_idx = [];
    for (let i = 0; i < p.length; i += 2) {
        let idx = add_point(p[i], p[i + 1]);
        p_idx.push(idx);
    }

    let p_2d_idx = p_idx.map(x => x * 2);
    let inv_side_lengths = [];
    let rest_lengths = [];
    let rest_lengths_sq = [];
    let edge_vec_idx = [];
    let normal_depths = [];

    // INIT
    for (let i = 0; i < p_2d_idx.length; i++) {
        let si = (i + 1 < p_2d_idx.length) ? i + 1 : 0; // successor of i

        // compute side lengths
        // TODO maybe just invert the given side_length arg
        let v = [
            points[p_2d_idx[si]] - points[p_2d_idx[i]],
            points[p_2d_idx[si] + 1] - points[p_2d_idx[i] + 1]
        ];
        inv_side_lengths[i] = 1.0 / Math.sqrt(v[0] * v[0] + v[1] * v[1]);

        // constraint length pairs
        edge_vec_idx.push(p_2d_idx[i], p_2d_idx[si]);
        let n = p_2d_idx.length - 3 - (i > 1 ? i - 1 : 0);
        for (let j = i + 2; j < i + n + 2; j++) {
            edge_vec_idx.push(p_2d_idx[i], p_2d_idx[j]);
        }

        // depth from each face
        let norm = [ // right-hand unit normal of i-th side
            (points[p_2d_idx[si] + 1] - points[p_2d_idx[i] + 1]) * inv_side_lengths[i],
            (points[p_2d_idx[i]]      - points[p_2d_idx[si]])    * inv_side_lengths[i]
        ];

        normal_depths[i] = 0;
        for (let j = 0; j < p_2d_idx.length; j++) {
            let proj = (points[p_2d_idx[j]] - points[p_2d_idx[i]]) * norm[0] + (points[p_2d_idx[j] + 1] - points[p_2d_idx[i] + 1]) * norm[1];
            normal_depths[i] = (normal_depths[i] < -proj) ? -proj : normal_depths[i];
        }
    }

    // constraint rest lengths
    for (let i = 0; i < edge_vec_idx.length / 2; i++) {
        let v = [
            points[edge_vec_idx[i * 2 + 1]]     - points[edge_vec_idx[i * 2]],
            points[edge_vec_idx[i * 2 + 1] + 1] - points[edge_vec_idx[i * 2] + 1]
        ];
        rest_lengths[i] = Math.sqrt(v[0] * v[0] + v[1] * v[1]);
        rest_lengths_sq[i] = rest_lengths[i] * rest_lengths[i];
    }

    return {
        p_idx,
        p_2d_idx,
        inv_side_lengths,
        rest_lengths,
        rest_lengths_sq,
        edge_vec_idx,
        normal_depths,
        point_generator,
        side_length,
        interacting: true,
        clipping_plane: null,
        hall_pass: null,
        color: rgba_to_hex(...color_map[num_sides]),
        fill_color: rgba_to_hex(...desaturate(...color_map[num_sides])),

        render() {
            context.beginPath();
            context.strokeStyle = '#' + this.color;
            context.lineWidth = 4;
            context.fillStyle = '#' + this.fill_color;
            for (let idx of this.p_2d_idx) {
                context.lineTo(points_transformed[idx], points_transformed[idx + 1]);
                //context.fillText(i, points_transformed[this.p_2d_idx[i]], canvas.height - points_transformed[this.p_2d_idx[i]+1]);
            }
            context.closePath();
            context.stroke();
            context.fill();
        },

        constrain() {
            for (let i = 0; i < this.edge_vec_idx.length / 2; i++) {
                let vecH = i * 2;
                let delta = [
                    points[this.edge_vec_idx[vecH + 1]]     - points[this.edge_vec_idx[vecH]],
                    points[this.edge_vec_idx[vecH + 1] + 1] - points[this.edge_vec_idx[vecH] + 1]
                ];
                let delSq = delta[0] * delta[0] + delta[1] * delta[1];
                let im1 = 1;//this.scene.pointsInvMass[this.edge_vec_idx[vecH]/2];
                let im2 = 1;//this.scene.pointsInvMass[this.edge_vec_idx[vecH+1]/2];
                //diff = (delSq-this.rest_lengths_sq[i])/((delSq+this.rest_lengths_sq[i])*(im1+im2)*2);
                let diff = (delSq - this.rest_lengths_sq[i]) / ((delSq + this.rest_lengths_sq[i]) * 4);
                for (let j = 0; j < 2; j++) {
                    let del = diff * delta[j];
                    points[this.edge_vec_idx[vecH] + j]   += del * im1;
                    points[this.edge_vec_idx[vecH+1] + j] -= del * im2;
                }
            }
        },

        collide_poly(poly2) {
            let minOverlap = 1000;
            let bestNorm = [];
            let besti = -1;
            let bestj = -1;
            let bestp = -1;

            // First check this poly's sides for penetration against "that" poly's points
            for (let i = 0; i < this.p_2d_idx.length; i++) {
                let j = (i + 1 < this.p_2d_idx.length) ? i + 1 : 0;
                let norm = [ // unit normal of i-th side
                    (points[this.p_2d_idx[j] + 1] - points[this.p_2d_idx[i] + 1]) * this.inv_side_lengths[i],
                    (points[this.p_2d_idx[i]]     - points[this.p_2d_idx[j]])     * this.inv_side_lengths[i]
                ];

                let min12Proj = 0;
                let max12Proj = 0;
                let min12Projp = -1;
                for (let p = 0; p < poly2.p_2d_idx.length; p++) {
                    let proj = (points[poly2.p_2d_idx[p]] - points[this.p_2d_idx[i]]) * norm[0] + (points[poly2.p_2d_idx[p] + 1] - points[this.p_2d_idx[i] + 1]) * norm[1];
                    if (proj < min12Proj) {
                        min12Proj = proj;
                        min12Projp = p;
                    }
                    max12Proj = (max12Proj < proj) ? proj : max12Proj;
                }

                if (min12Projp == -1) {
                    return false;
                }

                let minProj = (-this.normal_depths[i] < min12Proj) ? -this.normal_depths[i] : min12Proj;
                let maxProj = (max12Proj > 0) ? max12Proj : 0;

                let overlap = max12Proj - min12Proj + this.normal_depths[i] - maxProj + minProj;
                if (overlap > 0) {
                    if (overlap < minOverlap) {
                        minOverlap = overlap;
                        bestNorm = norm;
                        besti = this.p_2d_idx[i];
                        bestj = this.p_2d_idx[j];
                        bestp = poly2.p_2d_idx[min12Projp];
                    }
                } else {
                    return false;
                }
            }

            // Then check "that" poly's sides
            for (let i = 0; i < poly2.p_2d_idx.length; i++) {
                let j = (i + 1 < poly2.p_2d_idx.length) ? i + 1 : 0;
                let norm = [ // unit normal of i-th side
                    (points[poly2.p_2d_idx[j] + 1] - points[poly2.p_2d_idx[i] + 1]) * poly2.inv_side_lengths[i],
                    (points[poly2.p_2d_idx[i]]     - points[poly2.p_2d_idx[j]])     * poly2.inv_side_lengths[i]
                ];

                let min21Proj = 0;
                let max21Proj = 0;
                let min21Projp = -1;
                for (let p = 0; p < this.p_2d_idx.length; p++) {
                    let proj = (points[this.p_2d_idx[p]] - points[poly2.p_2d_idx[i]]) * norm[0] + (points[this.p_2d_idx[p] + 1] - points[poly2.p_2d_idx[i] + 1]) * norm[1];
                    if (proj < min21Proj) {
                        min21Proj = proj;
                        min21Projp = p;
                    }
                    max21Proj = (max21Proj < proj) ? proj : max21Proj;
                }

                if (min21Projp == -1) {
                    return false;
                }

                let minProj = (-poly2.normal_depths[i] < min21Proj) ? -poly2.normal_depths[i] : min21Proj;
                let maxProj = (max21Proj > 0) ? max21Proj : 0;

                let overlap = max21Proj - min21Proj + poly2.normal_depths[i] - maxProj + minProj;
                if (overlap > 0) {
                    //alert(overlap);
                    if (overlap < minOverlap) {
                        minOverlap = overlap;
                        bestNorm = norm;
                        besti = poly2.p_2d_idx[i];
                        bestj = poly2.p_2d_idx[j];
                        bestp = this.p_2d_idx[min21Projp];
                    }
                } else {
                    return false;
                }
            }

            // Update points
            let ijHitF;
            if (Math.abs(points[bestj] - points[besti]) > Math.abs(points[bestj + 1] - points[besti + 1])) {
                ijHitF = Math.abs(points[bestp] - points[besti]) / Math.abs(points[bestj] - points[besti]);
            } else {
                ijHitF = Math.abs(points[bestp + 1] - points[besti + 1]) / Math.abs(points[bestj + 1] - points[besti + 1]);
            }
            let l = 1.0 / (ijHitF * ijHitF + (1 - ijHitF) * (1 - ijHitF));
            let iHitP = (1 - ijHitF) * l;
            let jHitP = ijHitF * l;
            let push = [bestNorm[0] * minOverlap * 0.5, bestNorm[1] * minOverlap * 0.5];
            //alert (bestNorm[0] + ' ' + bestNorm[1]);
            points[besti]   -= push[0] * iHitP;
            points[besti+1] -= push[1] * iHitP;
            points[bestj]   -= push[0] * jHitP;
            points[bestj+1] -= push[1] * jHitP;
            points[bestp]   += push[0];
            points[bestp+1] += push[1];

            return true;
        },

        //PolyItem.prototype.snapPoly = function(poly2) {
        //    var points = this.scene.points;
        //    var pointsLast = this.scene.pointsLast;
        //    var pointsAccel = this.scene.pointsAccel;
        //
        //    // check all this poly's sides against all that poly's sides
        //    for (var i=0; i<this.p_2d_idx.length; i++) {
        //        var si = (i+1<this.p_2d_idx.length)?i+1:0;
        //
        //        for (var j=0; j<poly2.p_2d_idx.length; j++) {
        //            var sj = (j+1<poly2.p_2d_idx.length)?j+1:0;
        //
        //            var v1 = [points[poly2.p_2d_idx[sj]]-points[this.p_2d_idx[i]], points[poly2.p_2d_idx[sj]+1]-points[this.p_2d_idx[i]+1]];
        //            var v2 = [points[poly2.p_2d_idx[j]]-points[this.p_2d_idx[si]], points[poly2.p_2d_idx[j]+1]-points[this.p_2d_idx[si]+1]];
        //            var ds1 = v1[0]*v1[0]+v1[1]*v1[1]; // square of distance from this poly's point i to that poly's point sj
        //            var ds2 = v2[0]*v2[0]+v2[1]*v2[1];
        //
        //            if (ds1<100.0 && ds2<100.0) {
        //                var snapGap = 0.01;
        //                var ui = [(points[this.p_2d_idx[si]]-points[this.p_2d_idx[i]])*this.inv_side_lengths[i], (points[this.p_2d_idx[si]+1]-points[this.p_2d_idx[i]+1])*this.inv_side_lengths[i]];
        //                var uj = [(points[poly2.p_2d_idx[sj]]-points[poly2.p_2d_idx[j]])*poly2.inv_side_lengths[j], (points[poly2.p_2d_idx[sj]+1]-points[poly2.p_2d_idx[j]+1])*poly2.inv_side_lengths[j]];
        //                var ni = [ui[1], -ui[0]];
        //                var nj = [uj[1], -uj[0]];
        //                var v1 = [points[poly2.p_2d_idx[sj]]+nj[0]*snapGap-points[this.p_2d_idx[i]], points[poly2.p_2d_idx[sj]+1]+nj[1]*snapGap-points[this.p_2d_idx[i]+1]];
        //                var v2 = [points[poly2.p_2d_idx[j]]-points[this.p_2d_idx[si]]+nj[0]*snapGap, points[poly2.p_2d_idx[j]+1]-points[this.p_2d_idx[si]+1]+nj[1]*snapGap];
        //
        //                var vi = [points[this.p_2d_idx[i]]-pointsLast[this.p_2d_idx[i]], points[this.p_2d_idx[i]+1]-pointsLast[this.p_2d_idx[i]+1]];
        //                var vsi = [points[this.p_2d_idx[si]]-pointsLast[this.p_2d_idx[si]], points[this.p_2d_idx[si]+1]-pointsLast[this.p_2d_idx[si]+1]];
        //                var vj = [points[poly2.p_2d_idx[j]]-pointsLast[poly2.p_2d_idx[j]], points[poly2.p_2d_idx[j]+1]-pointsLast[poly2.p_2d_idx[j]+1]];
        //                var vsj = [points[poly2.p_2d_idx[sj]]-pointsLast[poly2.p_2d_idx[sj]], points[poly2.p_2d_idx[sj]+1]-pointsLast[poly2.p_2d_idx[sj]+1]];
        //
        //                pointsAccel[this.p_2d_idx[i]]     +=  0.0001*v1[0] - vi[0]*0.0001;
        //                pointsAccel[this.p_2d_idx[i]+1]   +=  0.0001*v1[1] - vi[1]*0.0001;
        //                pointsAccel[this.p_2d_idx[si]]    +=  0.0001*v2[0] - vsi[0]*0.0001;
        //                pointsAccel[this.p_2d_idx[si]+1]  +=  0.0001*v2[1] - vsi[1]*0.0001;
        //                pointsAccel[poly2.p_2d_idx[j]]    += -0.0001*v2[0] - vj[0]*0.0001;
        //                pointsAccel[poly2.p_2d_idx[j]+1]  += -0.0001*v2[1] - vj[1]*0.0001;
        //                pointsAccel[poly2.p_2d_idx[sj]]   += -0.0001*v1[0] - vsj[0]*0.0001;
        //                pointsAccel[poly2.p_2d_idx[sj]+1] += -0.0001*v1[1] - vsj[1]*0.0001;
        //
        //                return true;
        //            }
        //        }
        //    }
        //
        //    return false;
        //}

        point_intersects(point) {
            for (let i = 0; i < this.p_2d_idx.length; i++) {
                let j = (i + 1 < this.p_2d_idx.length) ? i + 1 : 0;
                let norm = [ // unit normal of i-th side
                    (points[this.p_2d_idx[j] + 1] - points[this.p_2d_idx[i] + 1]) * this.inv_side_lengths[i],
                    (points[this.p_2d_idx[i]]     - points[this.p_2d_idx[j]])     * this.inv_side_lengths[i]
                ];
                let proj = (point[0] - points[this.p_2d_idx[i]]) * norm[0] + (point[1] - points[this.p_2d_idx[i] + 1]) * norm[1];
                if (proj > 0) {
                    return false;
                }
            }

            return true;
        },

        // Checks the cross product of the first two sides. This assumes a convex poly.
        is_inverted() {
            let v1 = [
                points[this.p_2d_idx[1]]     - points[this.p_2d_idx[0]],
                points[this.p_2d_idx[1] + 1] - points[this.p_2d_idx[0] + 1]
            ];
            let v2 = [
                points[this.p_2d_idx[2]]     - points[this.p_2d_idx[1]],
                points[this.p_2d_idx[2] + 1] - points[this.p_2d_idx[1] + 1]
            ];
            return v1[0] * v2[1] - v1[1] * v2[0] < 0;
        },

        // Reverses the vertices' direction
        invert() {
            let inverted = [];
            for (let i = this.p_2d_idx.length - 1; i >= 0; i--) {
                inverted.push(this.p_2d_idx[i]);
            }
            this.p_2d_idx = inverted;
        },

        nudge_color() {
            let new_color = permute_color(...hex_to_rgba(this.color));
            this.color = rgba_to_hex(...new_color);
            this.fill_color = rgba_to_hex(...desaturate(...new_color));
        },

        clone() {
            let avg_x = this.p_idx.reduce((prev, curr) => prev + points[curr * 2], 0) / this.p_idx.length;
            let avg_y = this.p_idx.reduce((prev, curr) => prev + points[curr * 2 + 1], 0) / this.p_idx.length;
            let new_item = make_poly_item(this.point_generator, avg_x, avg_y, this.p_idx.length, this.side_length);
            new_item.nudge_color();
            new_item.clipping_plane = this.clipping_plane;
            poly_items.push(new_item);
            this.clipping_plane.add_poly_item(new_item);
            this.hall_pass = make_hall_pass(this.clipping_plane, 'clip_right', main_clipping_plane);

            return new_item;
        }
    };
}

function make_clipping_plane(clip_left, clip_right, clip_bottom, clip_top, hall_passes) {
    return {
        p_2d_idx: [],
        clip_left,
        clip_right,
        clip_bottom,
        clip_top,
        hall_passes,

        add_poly_item(poly_item) {
            for (let idx of poly_item.p_idx) {
                this.add_p_idx(idx);
            }
        },

        add_p_idx(...idx) {
            for (let i of idx) {
                this.p_2d_idx.push(i * 2);
            }
        },

        clip() {
            for (let idx of this.p_2d_idx) {
                if (points[idx] < this.clip_left) {
                    points[idx] = this.clip_left + Math.random() - 0.5;
                }
                // TODO somehow make a portal in the right clip plane for transfer to main clipping plane
                if (points[idx] > this.clip_right) {
                    points[idx] = this.clip_right + Math.random() - 0.5;
                }
                if (points[idx + 1] < this.clip_bottom) {
                    points[idx + 1] = this.clip_bottom + Math.random() - 0.5;
                }
                if (points[idx + 1] > this.clip_top) {
                    points[idx + 1] = this.clip_top + Math.random() - 0.5;
                }
            }
        },

        check_hall_passes() {
            
        }
    };
}

// TODO instead make portal passes consisting of p_idx and assign to clipping planes
function make_hall_pass(src, side, dest) {
    return {
        src,
        side,
        dest
    };
}

function render() {
    let half_width = canvas.width * 0.5;
    let half_height = canvas.height * 0.5;

    if (! matrix_good) {
        matrix.identity();
        matrix.scale(scale, scale);
        matrix.translate(half_width - pan_x, half_height - pan_y);
        //matrix.rotate(this.rotation);
        matrix_good = true;
    }

    points_transformed = matrix.transformArray(points);

    context.save();

    context.fillStyle = '#000000';
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.globalCompositeOperation = 'screen';

    for (let item of poly_items) {
        item.render();
    }

    palette.render();

    if (draw_text) {
        context.fillStyle = '#aaaaaa';
        context.textBaseline = 'top';
        context.font = '12pt Helvetica';
        context.fillText('frame rate: ' + frame_rate, 5, 5);
        context.fillText('collide rate: ' + collide_rate + ' / sec', 5, 20);
        context.fillText('polys: ' + poly_items.length, 5, 35);
    }

    context.restore();
}

function accumulate_forces() {
    for (let idx of palette.clipping_plane.p_2d_idx) {
        points_accel[idx + 1] = gravity_constant;
    }
}

function verlet(time_del) {
    for (let i = 0; i < num_points; i++) {
        let i2 = i * 2;
        let x_idx = i2;
        let y_idx = i2 + 1;
        let td2 = time_del * time_del;

        let temp_x = points[x_idx];
        let temp_y = points[y_idx];

        //points[x_idx] += (points[x_idx]-pointsLast[x_idx])*0.99 + pointsAccel[x_idx]*timeDel*timeDel;
        //points[y_idx] += (points[y_idx]-pointsLast[y_idx])*0.99 + pointsAccel[y_idx]*timeDel*timeDel;
        points[x_idx] += points[x_idx] - points_last[x_idx] + points_accel[x_idx] * td2;
        points[y_idx] += points[y_idx] - points_last[y_idx] + points_accel[y_idx] * td2;

        points_last[x_idx] = temp_x;
        points_last[y_idx] = temp_y;
    }
}

function constrain() {
    for (let item of poly_items) {
        item.constrain();
    }
    if (drag_item !== null) {
        drag_item.constrain();
    }
}

function collide_polys() {
    for (let i = 0; i < poly_items.length; i++) {
        for (let j = i + 1; j < poly_items.length; j++) {
            if (poly_items[i].interacting
                    && poly_items[j].interacting
                    && poly_items[i].clipping_plane === poly_items[j].clipping_plane) {

                let collision = poly_items[i].collide_poly(poly_items[j]);
                if (collision) {
                    collisions++;
                }
                //this.poly_items[i].color = '#ee0000';
                //this.poly_items[j].color = '#ee0000';
            }
        }
    }
}

function clip() {
    for (let plane of clipping_planes) {
        plane.clip();
    }
}

function safe_actions() {
    if (delete_drag) {
        if (drag_item !== null) {
            drag_item.poly_item.interacting = true;
            drag_item = null;
        }
        delete_drag = false;
    }

    for (let item of poly_items) {
        if (item.is_inverted()) {
            item.invert();
        }
    }
}

function tic() {
    render();
    accumulate_forces();
    verlet(render_interval);
    constrain();
    collide_polys();
    collide_polys();
    constrain();
    collide_polys();
    collide_polys();
    clip();
    safe_actions();
}

let color_map = {
    10: [255, 0, 0],//'#ff0000', // red
    9:  [255, 157, 0],//'#ff9d00', // orange
    8:  [255, 251, 0],//'#fffb00', // yellow
    7:  [128, 255, 0],//'#80ff00', // green
    6:  [0, 200, 255],//'#00c8ff', // cyan
    5:  [0, 64, 255],//'#0040ff', // blue
    4:  [128, 0, 255],//'#8000ff', // indigo
    3:  [221, 0, 255],//'#dd00ff'  // fuscia
}

function component_to_hex(c) {
    let hex = c.toString(16);
    return hex.length===1 ? '0' + hex : hex;
}

function rgba_to_hex(r, g, b, a = 255) {
    return component_to_hex(r) + component_to_hex(g) + component_to_hex(b) + component_to_hex(a);
}

function hex_to_rgba(hex) {
    let r = parseInt(hex.slice(0, 2), 16);
    let g = parseInt(hex.slice(2, 4), 16);
    let b = parseInt(hex.slice(4, 6), 16);
    let a = parseInt(hex.slice(6, 8), 16);
    return [r, g, b, a];
}

function bound(x, low = 0, high = 255) {
    return (x >= low) ? ((x <= high) ? x : high) : low;
}

function rgb_to_hsl(r, g, b) {
    let r_n = r / 255;
    let g_n = g / 255;
    let b_n = b / 255;

    let h, s, l;

    let c_max = Math.max(r_n, g_n, b_n);
    let c_min = Math.min(r_n, g_n, b_n);
    let del = c_max - c_min;

    l = (c_max + c_min) / 2;

    if (del === 0) {
        h = 0;
        s = 0;
    } else {
        s = del / (1 - Math.abs(2 * l - 1));

        if (c_max === r_n) {
            h = 60 * (((g_n - b_n) / del) % 6);
        } else if (c_max === g_n) {
            h = 60 * (((b_n - r_n) / del) + 2);
        } else {
            h = 60 * (((r_n - g_n) / del) + 4);
        }

        h = bound(h, 0, 359);
    }

    return [h, s, l];
}

function hsl_to_rgb(h, s, l) {
    let c = (1 - Math.abs(2 * l - 1)) * s;
    let x = c * (1 - Math.abs((h / 60) % 2 - 1));
    let m = l - c / 2;

    let r, g, b;

    switch (true) {
        case 0 <= h && h < 60:
            [r, g, b] = [c, x, 0];
            break;
        case 60 <= h && h < 120:
            [r, g, b] = [x, c, 0];
            break;
        case 120 <= h && h < 180:
            [r, g, b] = [0, c, x];
            break;
        case 180 <= h && h < 240:
            [r, g, b] = [0, x, c];
            break;
        case 240 <= h && h < 300:
            [r, g, b] = [x, 0, c];
            break;
        case 300 <= h && h < 360:
            [r, g, b] = [c, 0, x];
            break;
        default:
            console.log('Invalid hue: ' + String(h));
            [r, g, b] = [1, 1, 1];
            break;
    }

    return [r, g, b].map(x => Math.round((x + m) * 255));
}

function permute_color(r, g, b, a, range = 30) {
    let [h, s, l] = rgb_to_hsl(r, g, b);
    h = Math.round(h + Math.random() * range - range * 0.5);
    h = (h + 360) % 360;
    [r, g, b] = hsl_to_rgb(h, s, l);
    return [r, g, b];//.map(x => Math.round(bound(x + Math.random() * 60 - 30)));
}

function desaturate(r, g, b, factor = 0.1) {
    return [r, g, b].map(x => Math.round(x * factor));
}

function add_point(x, y) {
    points.push(x, y);
    points_last.push(x, y);
    points_accel.push(0, 0);
    return num_points++;
}

function make_regular_poly_coords(x, y, n, l = 20) {
    let t = 2 * Math.PI / n
    let sl = Math.sin(t / 2);
    let r = l / sl;

    let p = [];
    for (let i = 0; i < n; i++) {
        p[2 * i] = x + r * Math.cos(t / 2 + i * t);
        p[2 * i + 1] = y + r * Math.sin(t / 2 + i * t);
    }

    return p;
}

function add_regular_poly(x, y, n, clipping_plane = main_clipping_plane, pure_color = false) {

    let poly_item = make_poly_item(make_regular_poly_coords, x, y, n, 20);

    poly_item.clipping_plane = clipping_plane;

    if (! pure_color) {
        poly_item.nudge_color();
    }

    poly_items.push(poly_item);

    clipping_plane.add_poly_item(poly_item);
}

function add_rhomb_a(x, y, clipping_plane = main_clipping_plane) {
    let l = 100;
    let p1 = add_point(x + l * Math.cos(Math.PI / 10), y);
    let p2 = add_point(x, y + l * Math.sin(Math.PI / 10));
    let p3 = add_point(x - l * Math.cos(Math.PI / 10), y);
    let p4 = add_point(x, y - l * Math.sin(Math.PI / 10));

    clipping_plane.add_p_idx(p1, p2, p3, p4);

    let poly_item = make_poly_item([p1, p2, p3, p4]);
    poly_item.clipping_plane = clipping_plane;
    poly_items.push(poly_item);
}

function add_rhomb_b(x, y, clipping_plane = main_clipping_plane) {
    let l = 100;
    let p1 = add_point(x + l * Math.cos(Math.PI / 5), y);
    let p2 = add_point(x, y + l * Math.sin(Math.PI / 5));
    let p3 = add_point(x - l * Math.cos(Math.PI / 5), y);
    let p4 = add_point(x, y - l * Math.sin(Math.PI / 5));

    clipping_plane.add_p_idx(p1, p2, p3, p4);

    let poly_item = make_poly_item([p1, p2, p3, p4]);
    poly_item.clipping_plane = clipping_plane;
    poly_items.push(poly_item);
}

function make_drag_item(x, y, poly_item) {
    let rest_lengths_sq = [];
    let p_2d_idx = poly_item.p_2d_idx;
    for (let idx of p_2d_idx) {
        let v = [
            points[idx] - x,
            points[idx + 1] - y
        ];
        rest_lengths_sq.push(v[0] * v[0] + v[1] * v[1]);
    }
    return {
        x,
        y,
        p_2d_idx,
        poly_item,
        rest_lengths_sq,
        constrain() {
            for (let i = 0; i < this.p_2d_idx.length; i++) {
                let delta = [
                    points[this.p_2d_idx[i]] - this.x,
                    points[this.p_2d_idx[i] + 1] - this.y
                ];
                let del_sq = delta[0] * delta[0] + delta[1] * delta[1];
                let diff = (del_sq - this.rest_lengths_sq[i]) / (del_sq + this.rest_lengths_sq[i]);
                for (let j = 0; j < 2; j++) {
                    let del = diff * delta[j];
                    points[this.p_2d_idx[i] + j] -= del;
                }
            }
        }
    };
}

function make_palette_grip(rel_x) {
    return {
        rel_x
    };
}

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
              pan_x -= 5;
          }
          break;
        case 38:                                      //up arrow
          pan_y -= 5;
          break;
        case 39:                                      //right arrow
          if (event.shiftKey) {
          } else {
              pan_x += 5;
          }
          break;
        case 40:                                      //down arrow
          pan_y += 5;
          break;
        case 46:                                      //delete
          break;
        case 48:                                      //'0'
          add_regular_poly(px, py, 10)
          break;
        case 49: case 50:                             //'1'-'2'
          break;
        case 51: case 52: case 53: case 54: case 55: case 56: case 57: //'3'-'9'
          add_regular_poly(px, py, event.which - 48);
          break;
        case 61:                                      //'=' / '+'
          break;
        case 65:                                      //'a'
          add_rhomb_a();
          break;
        case 66:                                      //'b'
          add_rhomb_b();
          break;
        case 67:                                      //'c'
          break;
        case 68:                                      //'d'
          addDart(); //TODO
          break;
        case 69:                                      //'e'
          break;
        case 70:                                      //'f'
          add_regular_poly(px, py, Math.floor(Math.random() * 8 + 3));
          break;
        case 71:                                      //'g'
          gravity = !gravity;
          break;
        case 72:                                      //'h'
          break;
        case 73:                                      //'i'
          scale *= 1.1;
          matrix_good = false;
          break;
        case 76:                                      //'l'
          break;
        case 78:                                      //'n'
          draw_text = ! draw_text;
          break;
        case 79:                                      //'o'
          scale /= 1.1;
          matrix_good = false;
          break;
        case 80:                                      //'p'
          pause = !pause;
          break;
        case 82:                                      //'r'
          break;
        case 83:                                      //'s'
          break;
        case 85:                                      //'u'
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

function is_touch_device() {
    return 'ontouchstart' in window        // works on most browsers 
        || navigator.maxTouchPoints;       // works on IE10/11 and Surface
};

function transform_point(p) {
    let canvas_bound = canvas.getBoundingClientRect();
    return [
        -canvas_bound.width * 0.5 + p[0],
        -canvas_bound.height * 0.5 + p[1]
    ];
}

function set_point(p) {
    px = p[0];// - canvas_bound.left;
    py = p[1];// + canvas_bound.top;
}

function point_start(p) {
    set_point(p);

    if (palette.does_handle_intersect(p)) {
        if (palette.does_plus_intersect(p)) {
            palette.double_size();
        } else if (palette.does_minus_intersect(p)) {
            palette.halve_size();
        } else {
            let rel_p = palette.get_handle_point(p);
            palette.grip = make_palette_grip(rel_p[0]);
        }
    } else {
        for (let item of poly_items) {
            if (item.point_intersects(p)) {
                if (item.clipping_plane === palette.clipping_plane) {
                    let new_item = item.clone();
                    new_item.interacting = false;
                    drag_item = make_drag_item(p[0], p[1], new_item);
                    break;
                } else {
                    //item.interacting = false;
                    drag_item = make_drag_item(p[0], p[1], item);
                    break;
                }
            }
        }
    }
}

function point_move(p) {
    set_point(p);

    if (palette.grip !== null) {
        let width = p[0] - palette.grip.rel_x;
        palette.resize(width);
    } else if (drag_item !== null) {
        drag_item.x = p[0];
        drag_item.y = p[1];
    }
}

function point_end() {
    palette.grip = null;
    delete_drag = true;
}

if (false && is_touch_device()) {

    canvas.addEventListener('touchstart', function(event) {
        let point = transform_point([event.touches[0].clientX, event.touches[0].clientY]);
        point_start(point);
    });

    canvas.addEventListener('touchmove', function(event) {
        let point = transform_point([event.touches[0].clientX, event.touches[0].clientY]);
        point_move(point);
    });

    canvas.addEventListener('touchend', point_end);

} else {

    canvas.addEventListener('mousedown', function(event) {
        let point = transform_point([event.clientX, event.clientY]);
        point_start(point);
    });

    canvas.addEventListener('mousemove', function(event) {
        let point = transform_point([event.clientX, event.clientY]);
        point_move(point);
    });

    canvas.addEventListener('mouseup', point_end);
}

setInterval(function () {
    tic();
    frames++;
}, render_interval);

setInterval(function () {
    frame_rate = frames;
    collide_rate = collisions;
    frames = 0;
    collisions = 0;
}, 1000);

window.onresize = function() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
    main_clipping_plane.clip_left = -canvas.width * 0.5;
    main_clipping_plane.clip_right = canvas.width * 0.5;
    main_clipping_plane.clip_bottom = -canvas.height * 0.5;
    main_clipping_plane.clip_top = canvas.height * 0.5;
    palette.resize(palette.width);
    matrix_good = false;
}

let event = new Event('resize');
window.dispatchEvent(event);

let x = (palette.clipping_plane.clip_left + palette.clipping_plane.clip_right) * 0.5;
let y = (palette.clipping_plane.clip_bottom + palette.clipping_plane.clip_top) * 0.5;
for (let n of [3, 4, 5, 6, 7, 8, 9, 10]) {
    add_regular_poly(x, y, n, palette.clipping_plane, true);
}

//add_rhomb_a(palette.clipping_plane);
//add_rhomb_b(palette.clipping_plane);
