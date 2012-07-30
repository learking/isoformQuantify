from django.conf.urls import patterns, include, url

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    url('^getIsoformAbund', 'isoformQuantifyApp.getIsoformAbund.returnIsoformAbund'),
    url(r'^isoformQuantify/gene_(\w+)\.json$', 'isoformQuantifyApp.views.initializeGeneJSON'),
    url(r'^isoformQuantify/$', 'isoformQuantifyApp.views.index'),

    url(r'^isoformQuantify/testGeneLookUp/$', 'isoformQuantifyApp.views.testGeneLookUp'),
    # Examples:
    # url(r'^$', 'isoformQuantify.views.home', name='home'),
    # url(r'^isoformQuantify/', include('isoformQuantify.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
)
