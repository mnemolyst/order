var express = require('express');
var router = express.Router();
var requirejs = require('requirejs');

/* GET home page. */
router.get('/', function(req, res, next) {
    res.render('index', {
    });
});

module.exports = router;