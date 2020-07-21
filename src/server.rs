
use actix_files::NamedFile;

use actix_web::http::header::{ContentDisposition, DispositionType};
use actix_web::{HttpRequest, Result};
use std::path::PathBuf;
use clap::ArgMatches;

pub struct Vis {
    matches: &ArgMatches,
    range: StringRegion,
    mut list: Vec<(u64, Record)>,
    mut annotation: Vec<(u64, bed::Record)>,
    mut frequency: BTreeMap<u64, Vec<(u64, u32)>>,
    lambda: F,
}

async fn index(req: HttpRequest) -> Result<NamedFile> {
    let zoom: i64 = req.match_info().query("zoom").parse().unwrap();
    let path= req.match_info().query("filename");// .parse().unwrap();
    match NamedFile::open(format!("{}/{}.png", zoom, path)) {
        Ok(file) => Ok(file
            .use_last_modified(true)
            .set_content_disposition(ContentDisposition {
                disposition: DispositionType::Attachment,
                parameters: vec![],
            })),
        _ => {

            Ok(NamedFile::open(format!("{}/{}.png", zoom, path))?
                .use_last_modified(true)
                .set_content_disposition(ContentDisposition {
                    disposition: DispositionType::Attachment,
                    parameters: vec![],
                }))
        }
    }
}

#[actix_rt::main]
pub async fn server(matches: &ArgMatches, threads: u16) -> std::io::Result<()> {
    use actix_web::{web, App, HttpServer};

    HttpServer::new(|| App::new().route("/{zoom:.*}/{filename:.*}.png", web::get().to(index)))
        .bind("0.0.0.0:4000")?
        .workers(threads as usize)
        .run()
        .await
}
