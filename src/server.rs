
use actix_files::NamedFile;

use actix_web::http::header::{ContentDisposition, DispositionType};
use actix_web::{HttpRequest, Result};
use std::{sync::Mutex, path::PathBuf, cell::Cell};
use clap::ArgMatches;
use ghi::vis::bam_record_vis;
use crate::subcommands::bam_vis;



async fn index(req: HttpRequest) -> Result<NamedFile> {
    let zoom: i64 = req.match_info().query("zoom").parse().unwrap();
    let path: i64 = req.match_info().query("filename").parse().unwrap();// .parse().unwrap();
    match NamedFile::open(format!("{}/{}.png", zoom, path)) {
        Ok(file) => Ok(file
            .use_last_modified(true)
            .set_content_disposition(ContentDisposition {
                disposition: DispositionType::Attachment,
                parameters: vec![],
            })),
        _ => {
            let arg_vec = vec!["ghb", "vis", "-t", "1", "-r",  "parse"];
            // bam_record_vis(matches, string_range, list, ann, freq, |_| None)?;
            bam_vis(matches, 1);
            Ok(NamedFile::open(format!("{}/{}.png", zoom, path))?
                .use_last_modified(true)
                .set_content_disposition(ContentDisposition {
                    disposition: DispositionType::Attachment,
                    parameters: vec![],
                }))
        }
    }
}

/*
pub struct Vis<'a> {
    range: StringRegion,
    list: Mutex<Vec<(u64, Record)>>,
    annotation: Mutex<Vec<(u64, bed::Record)>>,
    index_list: Mutex<Vec<usize>>,
    suppl_list: Mutex<Vec<(&'a [u8], usize, usize, usize, usize)>>,
}
*/

#[actix_rt::main]
pub async fn server(app: &clap::App, matches: &ArgMatches, threads: u16) -> std::io::Result<()> {

    std::env::set_var("RUST_LOG", "actix_web=info");
    use actix_web::{web, App, HttpServer};

    // Create some global state prior to building the server
    #[allow(clippy::mutex_atomic)] // it's intentional.
    let counter1 = web::Data::new(Mutex::new(app));

    HttpServer::new(move || {
        let counter2 = Cell::new(app);
        App::new().data(counter2).route("/{zoom:.*}/{filename:.*}.png", web::get().to(index))
    }
    )
        .bind("0.0.0.0:4000")?
        .workers(threads as usize)
        .run()
        .await
}
