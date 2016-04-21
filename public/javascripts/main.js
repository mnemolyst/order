requirejs.config({
    baseUrl: 'public/javascripts/lib',
    paths: {
        app: '../app',
        jquery: 'jquery-2.2.0.min',
        //jqueryui: 'jquery-ui-1.11.4.custom/jquery-ui.min',
        //paper: 'paper-core',
    },
    shim: {
        'Matrix2D': {
            exports: 'Matrix2D',
        },
    },
});

requirejs(['order']);
